"""
collection of functions to lose data in xarray with different patterns (random, satellite-swath, ...)

    @author: Verena Bessenbacher
    @date: 13 06 2019
"""

import glob
import numpy as np
import random
import xarray as xr
import xesmf as xe
np.random.seed(0) # otherwise lazy evaluation has different missing patterns for lostmask and data

# DELETION
def listwise_deletion(data):
    data = np.where(np.any(np.isnan(data), axis=0), True, False) # wrap in xarray?
    return data

# MCAR
def set_random_fraction_to_nan_arr(data, frac_missing): 
    np.random.seed(0)
    #np.random.seed(0) # otherwise lazy evaluation has different missing patterns for lostmask and data
    data = data.copy() # SHOYER solved this for dask arrays
    data[np.random.rand(*data.shape) < frac_missing] = np.nan
    return data

def set_random_fraction_to_nan_bas(data, frac_missing):

    if isinstance(data, xr.core.dataset.Dataset):
        raise TypeError('function only defined for xarray DataArrays')

    data = xr.apply_ufunc(set_random_fraction_to_nan_arr, data, output_core_dims=[['latitude','longitude']], input_core_dims=[['latitude','longitude']], dask='parallelized', output_dtypes=[data.dtype], kwargs={'frac_missing':frac_missing})

    return data

def set_random_fraction_to_nan_veri(data, frac_missing):
    v,t,i,j = data.shape
    size = data.size
    k = int(frac_missing*size)
    np.random.seed(0)
    indices = random.sample(range(size),k)
    indices = np.unravel_index(indices, data.shape)
    data.values[indices[0], indices[1], indices[2], indices[3]] = np.nan

    return data

# MNAR
# patterns in monthly ESA-CCI-SM_combined monthly: vegetation missingness
# patterns in daily ESA-CCI-SM_combined daily: vegetation missingness and satellite swaths

def delete_sm_data(data, vod=None, rainforest=None, porosity=None, wetland_fraction=None, topographic_complexity=None, snowicedays=None):

    data = data.squeeze() # get rid of 1-length variable dimension
    data = data.rename({'longitude':'lon','latitude':'lat'}) # xesmf needs these coordinate names

    # extract masking in monthly ESA-CCI-SM_combined
    esaccismcombinedpath = '/net/exo/landclim/data/variable/soil-moisture/ESA-CCI-SM_combined/v04.4/0.25deg_lat-lon_1d/original/ancillary/'
    filenames = glob.glob(esaccismcombinedpath + '*.nc') # avail 1978 to 2017
    datamask = xr.open_mfdataset(filenames, combine='by_coords')
    datamask = datamask.to_array().load()

    #datamask = datamask.sm
    #datamask = ~np.isnan(datamask)
    ##datamask = datamask.rename({'lon':'longitude','lat':'latitude'})
    regridder = xe.Regridder(datamask, data, 'bilinear', reuse_weights=True)
    datamask = regridder(datamask)
    datamask = datamask.fillna(0) # marginally different land-sea mask and missing values on land lead to deleting points in data that are not above cf threshold. aim is if cf = 1 no data is deleted

    # mask data
    if vod is not None:
        data = data.where(datamask.loc['vod'] < vod)

    if rainforest is not None:
        data = data.where(datamask.loc['rainforest'] < rainforest)

    if porosity is not None: # not sure which one is bad
        raise AttributeError('check out if actually low or high porosity is bad')
        #data = data.where(datamask.loc['porosity'] < porosity)

    if wetland_fraction is not None:
        data = data.where(datamask.loc['wetland_fraction'] < wetland_fraction)

    if topographic_complexity is not None:
        data = data.where(datamask.loc['topographic_complexity'] < topographic_complexity)

    if snowicedays is not None:
        # this is in dask bec otherwise very slow IO
        esaccisnowicepath = '/net/exo/landclim/data/variable/soil-moisture/ESA-CCI-SM_combined/v04.4/0.25deg_lat-lon_1d/processed/netcdf/'
        filenames = glob.glob(esaccisnowicepath + '*-2???.nc')
        dataflag = xr.open_mfdataset(filenames, combine='by_coords')
        flag = dataflag.flag
        mask = flag.where(flag == 1, 0)
        monthmask = mask.resample(time='1m').sum(dim='time')
        monthmask = monthmask.resample(time='1MS').asfreq()
        regridder = xe.Regridder(monthmask, data, 'bilinear', reuse_weights=True)
        monthmask = regridder(monthmask)
        data = data.where(monthmask.load() < snowicedays)

    data = data.rename({'lon':'longitude','lat':'latitude'}) # rename data back
        
    return data

def delete_cloud_cover_deprecated(data, cloud_frac):

    if cloud_frac is not None:


        # extract cloud masking from ESA-CCI-Cloud cover
        import matplotlib.pyplot as plt
        esaccicloudcoverpath = '/net/exo/landclim/data/variable/temperature/ESA-CCI-Cloud/v2.0/0.5deg_lat-lon_1m/original/MODIS-TERRA/**/' # avail 2001 to 2014
        #testfile = '/net/exo/landclim/data/variable/temperature/ESA-CCI-Cloud/v2.0/0.5deg_lat-lon_1m/original/MODIS-TERRA/2000/200004-ESACCI-L3C_CLOUD-CLD_PRODUCTS-MODIS_TERRA-fv2.0.nc'
        #testfiles = '/net/exo/landclim/data/variable/temperature/ESA-CCI-Cloud/v2.0/0.5deg_lat-lon_1m/original/MODIS-TERRA/2000/20000[3-4]-ESACCI-L3C_CLOUD-CLD_PRODUCTS-MODIS_TERRA-fv2.0.nc'
        filenames = glob.glob(esaccicloudcoverpath + '*.nc', recursive=True)
        filenames = sorted(filenames)
        test = xr.open_mfdataset(filenames[:11] + filenames[13:], combine='nested', concat_dim='time')
        for t in range(test.cfc.shape[0]):
            print(t, np.isnan(test.cfc[t,:,:]).sum().values)
        quit()
        datamask = xr.open_mfdataset(filenames, combine='by_coords')
        cfc = datamask.cfc
        del(datamask)
        datamask = datamask.where(~np.isnan(datamask),1.)
        datamask = datamask < cloud_frac
        datamask = datamask.resample(time='1MS').asfreq().astype(bool) # sync time stamp

        data = data.rename({'longitude':'lon','latitude':'lat'}) # xesmf needs these coordinate names
        regridder = xe.Regridder(datamask, data, 'bilinear', reuse_weights=True)
        datamask = regridder(datamask)

        data = data.where(datamask)
        data = data.rename({'lon':'longitude','lat':'latitude'}) # rename data back

    return data

def delete_cloud_cover(data, cloud_frac):

    if cloud_frac is not None:

        data = data.rename({'longitude':'lon','latitude':'lat'}) # xesmf needs these coordinate names

        # extract cloud masking from ESA-CCI-Cloud cover
        esaccicloudcoverpath = '/net/exo/landclim/data/variable/temperature/ESA-CCI-Cloud/v2.0/0.5deg_lat-lon_1m/original/MODIS-TERRA/**/' # avail 2001 to 2014
        filenames = glob.glob(esaccicloudcoverpath + '*.nc', recursive=True)
        filenames = sorted(filenames)
        filenames = filenames[11:] # excluding year 2000

        datamask_all = np.full_like(data, np.nan)

        for f, filename in enumerate(filenames):
            cfc = xr.open_dataset(filename).cfc
            datamask = cfc < cloud_frac
            datamask = datamask.resample(time='1MS').asfreq().astype(bool) # sync time stamp

            regridder = xe.Regridder(datamask, data, 'bilinear', reuse_weights=True)
            datamask_all[f,:,:] = regridder(datamask)

        data = data.where(datamask_all)
        data = data.rename({'lon':'longitude','lat':'latitude'}) # rename data back

    return data
    

def delete_rain_mask(data, imask=True, relweight=None):

    precipfiles = '/home/bverena/local_data_store/GPM/*.nc'
    filenames = glob.glob(precipfiles, recursive=True)
    precipdata = xr.open_mfdataset(filenames, combine='by_coords')
    precipdata['time'] = precipdata.time.to_index().to_datetimeindex(unsafe=True)
    precipdata = precipdata.resample(time='1MS').asfreq() # sync time stamp
    # this does not work since CFTimeIndex cannot be converted silently to datetimeindex. There is a method xarray.CFTimeIndex.to_datetimeindex but I cannot find the xr.CFTimeIndex object within precipdata. possibly related issue fix https://github.com/pydata/xarray/issues/2191 is too slow causing python internal error.; update: found the object for the method. it works and warning can be ignored because I checked the times are correct
    datamask = precipdata['precipitation'].load()
    # from documentation:
    # https://pmm.nasa.gov/sites/default/files/document_files/IMERG_doc_190313.pdf
    # "the  interpretation  is  that  this  is  the  approximate  number  of  gauges  required  to  produce  the estimated random error, given the estimated precipitation.estimated random error, given the estimated precipitation"
    # my interpretation: if the random error is small, we need many gauges to produce a similarly small random error. i.e. if gaugeRelativeWeighting is large, the error is small.
    # no it is exactly the other way around 13 12 2019
    #dataqf = precipdata['gaugeRelativeWeighting'].load()
    dataqf = precipdata['precipitationQualityIndex'].load()

    datamask = datamask.transpose('time', 'lat', 'lon')
    dataqf = dataqf.transpose('time', 'lat', 'lon')

    data = data.rename({'longitude':'lon','latitude':'lat'}) # xesmf needs these coordinate names
    regridder = xe.Regridder(datamask, data, 'bilinear', reuse_weights=True) #TODO numpy.ascontiguousarray(x)
    datamask = regridder(datamask)
    regridder = xe.Regridder(dataqf, data, 'bilinear', reuse_weights=True)
    dataqf = regridder(dataqf)

    if imask:
        #data.loc[datamask.time,:,:] = data.where(datamask)
        data = data.where(~np.isnan(datamask))
    if relweight is not None:
        #data.loc[datamask.time,:,:] = data.where(dataqf > relweight*100)
        #data = data.where(~np.isnan(dataqf))
        #data = data.where(~np.isnan(dataqf)) # this does not do anything currently
        data = data.where(dataqf > relweight)
    data = data.rename({'lon':'longitude','lat':'latitude'}) # rename data back

    return data


if __name__ == "__main__":

    # specify file path
    era5path = '/net/exo/landclim/data/dataset/ERA5_deterministic/recent/0.25deg_lat-lon_1m/original/'
    filenames = glob.glob(era5path + 'era5_deterministic_recent.*.025deg.1m.2014.nc')

    # filter for files with variable names in varlist
    varlist = ['ro', 'tp', 'e', 'swvl1', 't2m', 'slt', 'lai_hv', 'swvl2', 'swvl3', 'swvl4', 'skt'] 
    filenames = [filename for filename in filenames if any(varname == filename.split('.')[2] for varname in varlist)]

    # load file
    dataxr = xr.open_mfdataset(filenames, combine='by_coords')
    dataxr = dataxr.to_array().load()

    # set fraction missing
    frac_missing = 0.5

    # lose data MNAR according to vegetation density
    print('start deleting tp ...')
    delete_rain_mask(dataxr.loc['tp'])
    dataxr.loc['tp'] = delete_rain_mask(dataxr.loc['tp'])
    print('start deleting cc ...')
    dataxr.loc['skt'] = delete_cloud_cover(dataxr.loc['skt'], 0.5)
    print('start deleting sm ...')
    dataxr.loc['swvl1'] = delete_sm_data(dataxr.loc['swvl1'])
    #import matplotlib.pyplot as plt
    #cmap = plt.cm.PuOr
    #cmap.set_bad('grey')
    #for i in range(11):
    #    dataxr[i,0,:,:].plot(cmap=cmap)
    #    plt.show()
    
    #import matplotlib.pyplot as plt
    #dataxr[0,:,:].plot()
    #plt.show()
    quit()

    # lose 50% of the data at random
    set_random_fraction_to_nan_veri(dataxr.to_array(), frac_missing)
    # %timeit: 9.55 s ± 86.9 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    xr.apply_ufunc(set_random_fraction_to_nan_arr, dataxr, output_core_dims=[['latitude','longitude']], input_core_dims=[['latitude','longitude']], dask='parallelized', output_dtypes=[dataxr.to_array().dtype], kwargs={'frac_missing':frac_missing})
    # %timeit: 775 µs ± 1.28 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
    set_random_fraction_to_nan_bas(dataxr.to_array(), frac_missing)
    # %timeit: 893 µs ± 2.84 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
