"""
ADD DESCRIPTION

@author: Verena Bessenbacher
@date: 13 05 2019
"""

# imports
import numpy as np
import numpy.ma as ma
import xarray as xr
import xarray.ufuncs as xu
import glob
import logging
logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.DEBUG)
logging.getLogger('matplotlib').setLevel(logging.WARNING) # no matplotlib logging
logging.getLogger('regionmask').setLevel(logging.WARNING) # no matplotlib logging
logging.getLogger('fiona').setLevel(logging.WARNING) # no matplotlib logging
import regionmask
import warnings
import matplotlib.pyplot as plt
from datetime import datetime
from imputethis import Imputation
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from anchorRegression import anchorRegression

# own functions import
from lose_data import set_random_fraction_to_nan_bas as lose_data
from lose_data import delete_sm_data, delete_cloud_cover, delete_rain_mask
from era_functions import calc_frac_missing, mean_impute, locf_impute, temporal_interp, spatial_mean

# imputation settings
plotpath = '/home/bverena/scripts_metalina/'
cf =  0.09 # larger cf --> less missing values !
maxiter = 10000
missingness = 'realistic_fast'
fill_method = 'anchorRegression'

# operational settings
idebugspace = True # only every 10th latlon point for faster computation

# methodological settings
ilogscaleprecip = False #only set True again if nans do not disappear and zero precip problem is solved
inormalise = True
ispatiotemporal = True
ideseasonalise = True
inoantarctica = True
inogreenland = True
ilatlondata = True
itimedata = True

# plot settings
isavenc = True
isave = True

# parameter settings for imputation methods
# only works for svdimpute right now
k = None # for SVDimpute
alpha = 1000 # for SVDimpute and RidgeRegression
epsilon = 1e-5 # for all # 1e-5 == 0.00001

# following all for missforest
n_trees = 50
max_depth = 10
min_samples_leaf = 10
min_samples_split = 10
max_features = 'auto'
bootstrap = 'True'
min_samples_leaf = 10000 #faster
min_samples_split = 10000 #faster
n_trees = 50 #faster

# timeit
start = datetime.now()

# settings for realistic missingness pattern
imask = False; relweight = (1- cf)*10.
vod = cf ; rainforest = cf ; wetland_fraction=cf; topographic_complexity = cf*100.; snowicedays = cf*20
cloud_fraction = cf

# list of variables to use
invarnames = ['lsm','z','slor','cvl','cvh', 'tvl', 'tvh']
varnames = ['tp', 'skt', 'swvl1']
#import IPython; IPython.embed()
# close with writing quit()

# filepaths for ERA5
era5path_monthly = '/net/exo/landclim/data/dataset/ERA5_deterministic/recent/0.25deg_lat-lon_1m/original/'
era5path_invariant = '/net/exo/landclim/data/dataset/ERA5_deterministic/recent/0.25deg_lat-lon_time-invariant/original/'

# catch filenames
years = list(np.arange(2001,2015))
filenames_var = ['{}era5_deterministic_recent.{}.025deg.1m.{}.nc'.format(era5path_monthly,varname,year) for year in years for varname in varnames]
filenames_invar = ['{}era5_deterministic_recent.{}.025deg.time-invariant.nc'.format(era5path_invariant,varname) for varname in invarnames]

# open files
logging.info('open files ...')
dataxr = xr.open_mfdataset(filenames_var, combine='by_coords').to_array().load() # no dask bec sucks
timeinvariant = xr.open_mfdataset(filenames_invar, combine='by_coords').to_array().load()

# extract land-sea mask
logging.info('extract land sea mask ...')
lsm = timeinvariant.loc['lsm']
invarnames.remove('lsm')
timeinvariant = timeinvariant.loc[invarnames] #Dataarray method
landmask = (lsm.squeeze() > 0.8).load() # land is 1, ocean is 0

# convert z to topo by dividing with g
# description here https://confluence.ecmwf.int/display/CKB/ERA5%3A+surface+elevation+and+orography
logging.info('convert z to topo ...')
timeinvariant.loc['z'] = timeinvariant.loc['z'] / 9.80665 

# remove Antarctica because of known soil moisture reanalysis problems there (source: Bas)
if inoantarctica:
    dataxr = dataxr.sel(latitude=slice(90,-60))
    timeinvariant = timeinvariant.sel(latitude=slice(90,-60))
    landmask = landmask.sel(latitude=slice(90,-60))

# smaller dataset for debugging
if idebugspace:
    logging.info('shrink data for testing ...')
    ith = 10
    #ith = 50 # for param fitting
    dataxr = dataxr[:,:,::ith,::ith]
    timeinvariant = timeinvariant[:,:,::ith,::ith]
    landmask = landmask[::ith,::ith]

if inogreenland:
    n_greenland = regionmask.defined_regions.natural_earth.countries_110.map_keys('Greenland')
    mask = regionmask.defined_regions.natural_earth.countries_110.mask(dataxr, lon_name='longitude', lat_name='latitude', wrap_lon=True)
    landmask = np.where(mask != n_greenland, landmask, False)
    icemask = (mask == n_greenland)

landlat, landlon = np.where(landmask)
frac_ocean = np.isnan(dataxr).sum().values / dataxr.size
shape = dataxr.shape

# lose data
logging.info(f'lose data with missingness {missingness} ...')
if missingness == 'random':
    cfs = [0.99, 0.69, 0.59, 0.09]
    frac_mis = [0.08, 0.31, 0.52, 0.72]
    dataxr_lost = lose_data(dataxr.copy(deep=True), frac_mis[cfs == cf])
elif missingness == 'realistic':

    dataxr_lost = dataxr.copy(deep=True)

    dataxr_lost.loc['skt'] = delete_cloud_cover(dataxr_lost.loc['skt'], cloud_fraction) 
    dataxr_lost.loc['swvl1'] = delete_sm_data(dataxr_lost.loc['swvl1'], vod=vod, 
                                              rainforest=rainforest, wetland_fraction=wetland_fraction, 
                                              topographic_complexity=topographic_complexity, 
                                              snowicedays=snowicedays) # only first layer
    dataxr_lost.loc['tp'] = delete_rain_mask(dataxr_lost.loc['tp'], imask=imask, relweight=relweight)

    if isave:
        np.isnan(dataxr_lost).to_netcdf(plotpath + 'pattern_missing_{}_{}.nc'.format(cf, idebugspace))
elif missingness == 'realistic_fast':
    mask = xr.open_dataarray(plotpath + 'pattern_missing_{}_{}.nc'.format(cf, idebugspace))
    dataxr_lost = dataxr.where(~mask)
else:
    raise AttributeError('missingness pattern {} not recognized'.format(missingness))
frac_missing = calc_frac_missing(dataxr_lost, frac_ocean)

# initial gapfill
logging.info(f'impute data with init gapfill ...')
if fill_method == 'listwise':
    dataxr_filled = dataxr_lost.where(~np.any(np.isnan(dataxr_lost)))
    timeinvariant = timeinvariant.where(~np.any(np.isnan(timeinvariant)))
elif fill_method == 'mean_impute':
    dataxr_filled = mean_impute(dataxr_lost)
elif fill_method == 'locf':
    dataxr_filled = locf_impute(dataxr_lost)
elif fill_method in ['missforestmean', 'missforestlocf', 'ridgemean', 
                     'ridgelocf', 'svdimputemean','svdimputelocf','gp', 'anchorRegression']:

    logging.info(f'{fill_method}: initial gapfill ...')
    if fill_method in ['missforestmean', 'ridgemean', 'svdimputemean', 'anchorRegression']:
        dataxr_filled = mean_impute(dataxr_lost)
    elif fill_method in ['missforestlocf','ridgelocf', 'svdimputelocf','gp']:
        dataxr_filled = locf_impute(dataxr_lost)
    else:
        raise AttributeError(f'fill_method {fill_method} not recognized')
else:
    raise AttributeError(f'fill_method {fill_method} not recognized')

# include spatial neighborhood in data
if ispatiotemporal:
    logging.info('include spatial & temporal neighborhood ...')
    temporal = xr.full_like(dataxr_filled,np.nan)
    spatial = xr.full_like(dataxr_filled,np.nan)
    for varname in varnames:
        temporal.loc[varname] = temporal_interp(dataxr_filled.sel(variable=varname))
        spatial.loc[varname] = spatial_mean(dataxr_filled.sel(variable=varname))
    temporal = temporal.to_dataset(dim='variable').rename({'skt': 'skt_t', 'swvl1': 'swvl1_t', 'tp':'tp_t'}).to_array()
    spatial = spatial.to_dataset(dim='variable').rename({'skt': 'skt_sp', 'swvl1': 'swvl1_sp', 'tp':'tp_sp'}).to_array()
    dataxr_lost = xr.concat([dataxr_lost, temporal, spatial], dim='variable')
    dataxr_filled = xr.concat([dataxr_filled, temporal, spatial], dim='variable')

# add latitude and longitude as predictors
if ilatlondata:
    logging.info('add latitude and longitude as predictors ...')
    tmp = timeinvariant.to_dataset(dim='variable')
    londata, latdata = np.meshgrid(timeinvariant.longitude, timeinvariant.latitude)
    tmp['latdata'] = (('latitude', 'longitude'), latdata)
    tmp['londata'] = (('latitude', 'longitude'), londata)
    timeinvariant = tmp.to_array()

# remove ocean points 
logging.info('remove ocean points ...')
dataxr = dataxr.isel(longitude=xr.DataArray(landlon, dims='landpoints'), 
                     latitude=xr.DataArray(landlat, dims='landpoints'))
dataxr_lost = dataxr_lost.isel(longitude=xr.DataArray(landlon, dims='landpoints'), 
                                        latitude=xr.DataArray(landlat, dims='landpoints'))
dataxr_filled = dataxr_filled.isel(longitude=xr.DataArray(landlon, dims='landpoints'), 
                                        latitude=xr.DataArray(landlat, dims='landpoints'))
timeinvariant = timeinvariant.isel(longitude=xr.DataArray(landlon, dims='landpoints'), 
                                        latitude=xr.DataArray(landlat, dims='landpoints'))

# remove seasonality
if ideseasonalise:
    logging.info('deseasonalise data ...')
    # new: seasonality from perfect data
    seasonality = dataxr.groupby('time.month').mean(dim='time')
    dataxr = dataxr.groupby('time.month') - seasonality
    dataxr_lost = dataxr_lost.groupby('time.month') - seasonality
    dataxr_filled = dataxr_filled.groupby('time.month') - seasonality


# log scale on precipitation
if ilogscaleprecip:
    logging.info('log scale precip ...')
    raise AttributeError('not yet implemented')

# add time as predictor
if itimedata:
    logging.info('add time as predictor ...')
    _, ntimesteps, nlandpoints = dataxr_filled.shape
    tmp = dataxr_filled.to_dataset(dim='variable')
    timedat = np.arange(ntimesteps)
    timedat = np.tile(timedat, nlandpoints).reshape(nlandpoints,*timedat.shape).T
    tmp['timedat'] = (('time','landpoints'), timedat)
    dataxr_filled = tmp.to_array()

if inormalise: #mean zero unit (1) standard deviation
    logging.info('normalise data ...')
    datamean = dataxr_filled.mean(dim=('time','landpoints'))
    datastd = dataxr_filled.std(dim=('time','landpoints'))
    dataxr_lost = (dataxr_lost - datamean) / datastd
    dataxr = (dataxr - datamean) / datastd
    dataxr_filled = (dataxr_filled - datamean) / datastd
    invarmean = timeinvariant.mean(dim=('time','landpoints'))
    invarstd = timeinvariant.std(dim=('time','landpoints'))
    timeinvariant = (timeinvariant - invarmean) / invarstd

timebeforefill = datetime.now()

# obtain lostmask
lostmask = xu.isnan(dataxr_lost) # missing (lost + orig missing + ocean) is True, rest is False

# impute data
logging.info(f'impute data with method {fill_method} ...')
if fill_method in ['mean_impute','locf']:
    pass
elif fill_method in ['missforestmean', 'missforestlocf', 'ridgemean', 
                     'ridgelocf', 'svdimputemean','svdimputelocf','gp','anchorRegression']:

    #  stack var and invar
    logging.info(f'{fill_method}: stack ...')
    ntimesteps = dataxr_filled.coords['time'].size
    timeinvariant = np.repeat(timeinvariant, ntimesteps, axis=1)
    timeinvariant['time'] = dataxr_filled['time']
    dataall_filled = xr.concat([dataxr_filled, timeinvariant], dim='variable')

    # flatten time and space points
    dataall_filled = dataall_filled.stack(datapoints=('time','landpoints'))
    lostmask = lostmask.stack(datapoints=('time','landpoints'))

    # impute missing data
    logging.info(f'{fill_method}: impute values iteratively ...')
    if fill_method in ['missforestmean','missforestlocf']:
        impute = Imputation(epsilon, RandomForestRegressor, dimension_reduction=None, maxiter=maxiter,
                                      max_depth=max_depth, n_estimators=n_trees, 
                                      min_samples_leaf=min_samples_leaf, 
                                      min_samples_split=min_samples_split, max_features=max_features, 
                                      bootstrap=bootstrap, n_jobs=15, verbose=0, oob_score=True)
    elif fill_method in ['ridgemean', 'ridgelocf']:
        impute = Imputation(epsilon, Ridge, dimension_reduction=None, alpha=alpha, maxiter=maxiter)
    elif fill_method in ['svdimputemean', 'svdimputelocf']:
        impute = Imputation(epsilon, Ridge, PCA, k=None, alpha=alpha, maxiter=maxiter)
    elif fill_method in ["anchorRegression"]:
        impute = Imputation(1e-4, anchorRegression, dimension_reduction=None, tradeoff_param = 1, maxiter=100) 
    if ispatiotemporal:
        impute.set_filter(spatial_mean, temporal_interp, (landlat, landlon), shape)

    dataxr_filled = impute.impute(dataall_filled.T, lostmask.T, varnames)

    # unflatten == reshape
    logging.info(f'{fill_method}: unstack ...')
    dataxr_filled = dataxr_filled.unstack('datapoints')
    lostmask = lostmask.unstack('datapoints')

else:
    raise AttributeError('fill method {} not recognized'.format(fill_method))

# drop constant variables
dataxr_filled = dataxr_filled.sel(variable=varnames)
dataxr_lost = dataxr_lost.sel(variable=varnames)
dataxr = dataxr.sel(variable=varnames)
lostmask = lostmask.sel(variable=varnames)
datamean = datamean.sel(variable=varnames)
datastd = datastd.sel(variable=varnames)

# calculate expensiveness
logging.info('calculate metrics ...')
timeafterfill = datetime.now()
timeit = timeafterfill - timebeforefill

# renormalise, exp precip
if inormalise:
    dataxr_filled = dataxr_filled * datastd + datamean
    dataxr_lost = dataxr_lost * datastd + datamean
    dataxr = dataxr * datastd + datamean

if fill_method in ['svdimputelocf', 'missforestlocf', 'ridgelocf', 
                   'svdimputemean', 'missforestmean', 'ridgemean', 'anchorRegression']:

    dataxr['landpoints'] = dataxr_filled['landpoints'] # in dataxr still with lat lon coords #only svdimpute?
    lostmask['landpoints'] = dataxr_filled['landpoints'] # in dataxr still with lat lon coords

# print and save metrics
print('********************* CALCULATION FINISHED *********************')
print('frac_missing: ', frac_missing, ' imputation method: ', fill_method, 'missingness pattern: ', missingness)
print('cf: ', cf, ' idebugspace: ', idebugspace, 'ispatiotemporal: ', ispatiotemporal)
print('timeit: ', timeit)

if isavenc:
    dataxr_filled.to_netcdf(plotpath + f'dataxr_filled_{missingness}_{cf}_{fill_method}_idebug_{idebugspace}_ispatiot_{ispatiotemporal}.nc')
    dataxr_lost.to_netcdf(plotpath + f'dataxr_lost_{missingness}_{cf}_idebug_{idebugspace}.nc')
    dataxr.to_netcdf(plotpath + f'dataxr_idebug_{idebugspace}.nc')
    np.save(plotpath + f'landmask_idebug_{idebugspace}', landmask)
    np.save(plotpath + f'icemask_idebug_{idebugspace}', icemask)
