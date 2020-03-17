"""
Class for iterative Imputation of missing values in geoscientific datasets. Datasets must be xarray with the shape (variables, timesteps, latitude, longitude), but flattened to (variables, datapoints). Returns same feature table but with estimation of missing values filled in. Cannot deal with missing values itself, therefore the missing values need to be gapfilled beforehand and the missingness needs to be passed via a boolean mask of the same shape as the data. DOCUMENTATION TBC

    @author: Verena Bessenbacher
    @date: 29 11 2019
"""

from numpy import Inf
import numpy as np
import logging
import warnings
logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.DEBUG)
logging.getLogger('matplotlib').setLevel(logging.WARNING) # no matplotlib logging
from anchorRegression import anchorRegression
import math
import random
import matplotlib.pyplot as plt
import copy
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from itertools import compress
import sys
class Imputation:
    """
    DOCSTRING
    """

    def __init__(self, epsilon, regression, dimension_reduction=None, maxiter=Inf, k=None, **kwargs):
        """
        regression: Class object that is a regression function, 
                    e.g. sklearn.ensemble.RandomForestRegressor, sklearn.linear_model.Ridge
        dimension_reduction: [optional] Class object that is any unsupervised machine learning method,
                    i.e. a dimension reduction, e.g. sklearn.decomposition.PCA
        k: no of dimensions retained after PCA. None equals all
        epsilon: convergence threshold for iteration
        **kwargs: additional arguments to be passed to regression initialisation
        """

        self.sup = regression(**kwargs)
        if dimension_reduction is not None:
            self.unsup = dimension_reduction(n_components=k)
        else:
            self.unsup = None
        self.epsilon = epsilon
        self.ispatiotemporal = False
        self.maxiter = maxiter

    def set_filter(self, spatial_filter, temporal_filter, idxs, shape):
        """
        spatial_filter: function like scipy.ndimage.filters.generic_filter with footprint in spatial dims
        temporal_filter: --"-- in temporal dims
        idxs: tuple of latitude and longitude indizes of original 4D data cube
        shape: tuple describing shape of original dataset without variables (ntime, nlat, nlon)
        """

        self.ispatiotemporal = True
        self.landlat, self.landlon = idxs
        self.shape = shape[1:]
        self.spatial_filter = spatial_filter
        self.temporal_filter = temporal_filter

    def _update_neighborhood_par(self, data, varnames):
        """
        data: feature table (stacked xarray) that can be unstacked along landpoints
        varnames: list of strings containing variable names where interpolation should take place
        """

        if len(varnames) != 3:
            raise AttributeError('_update_neighborhood_par function not yet implemented for other than 3 variables')

        # inflate dataset back to 4D to infer neighborhood
        tmp = np.full((data.shape[1],*self.shape), np.nan)
        tmp[:,:,self.landlat,self.landlon] = data.unstack()

        # TODO this is embarassingly parallel
        from concurrent.futures import ProcessPoolExecutor
        pool = ProcessPoolExecutor(6)

        n_v0 = np.where(data.coords['variable'] == varnames[0])[0][0] # TODO more elegant
        n_v1 = np.where(data.coords['variable'] == varnames[1])[0][0]
        n_v2 = np.where(data.coords['variable'] == varnames[2])[0][0]

        with warnings.catch_warnings():
            future0 = pool.submit(self.spatial_filter, tmp[n_v0,:,:])
            future1 = pool.submit(self.spatial_filter, tmp[n_v1,:,:])
            future2 = pool.submit(self.spatial_filter, tmp[n_v2,:,:])

            future3 = pool.submit(self.temporal_filter, tmp[n_v0,:,:])
            future4 = pool.submit(self.temporal_filter, tmp[n_v1,:,:])
            future5 = pool.submit(self.temporal_filter, tmp[n_v2,:,:])

        data.loc[:,f'{varnames[0]}_sp'] = future0.result()[:,self.landlat,self.landlon].flatten()
        data.loc[:,f'{varnames[1]}_sp'] = future1.result()[:,self.landlat,self.landlon].flatten()
        data.loc[:,f'{varnames[2]}_sp'] = future2.result()[:,self.landlat,self.landlon].flatten()

        data.loc[:,f'{varnames[0]}_t'] = future3.result()[:,self.landlat,self.landlon].flatten()
        data.loc[:,f'{varnames[1]}_t'] = future4.result()[:,self.landlat,self.landlon].flatten()
        data.loc[:,f'{varnames[2]}_t'] = future5.result()[:,self.landlat,self.landlon].flatten()

        # overwrite "island" nans with current best estimate
        for varname in varnames:
            data.loc[:,f'{varname}_sp'][np.isnan(data.loc[:,f'{varname}_sp'])] = data.loc[:,f'{varname}'][np.isnan(data.loc[:,f'{varname}_sp'])]

        return data

    def impute(self, data, lostmask, varnames, iremove=False):
        """
        data: feature table (stacked xarray) that can be unstacked along datapoints
        lostmask: xarray in the same shape as data, indicating which values were originally missing as True
        varnames: list of variables where estimates need to be updated
        iremove: remove all pieces of code that assume more about the data than given by vars. necessary to run example code in __main__ below
        """

        # define convergence criteria
        self.delta_n = Inf
        data_before = data.copy(deep=True) #DEBUG
        iverbose = False #DEBUG

        # set parameter settings
        self.n_miniter = 15 # minimal number of total iterations
        self.n_miniter_below = 5 # minimal number of iterations where epsilon is below threshold
        self.iter = 0
        self.iter_below = 0

        # while convergence not reached, loop over features and impute iteratively
        while True:

            # store previously imputed matrix
            data_old = data.copy(deep=True)

            for varname in varnames:

                # divide into predictor and predictands
                if iverbose:
                    logging.info('divide into predictor and predictands')
                y = data.loc[:,varname]
                notyvars = data.coords['variable'].values.tolist()
                notyvars.remove(varname)
                X = data.loc[:,notyvars]
                y_mask = lostmask.loc[:,varname]

                # fit dimension_reduction
                if self.unsup is not None:
                    X = self.unsup.fit_transform(X)

                # divide into missing and not missing values
                if iverbose:
                    logging.info('divide into missing and not missing')
                y_mis = y[y_mask]
                y_obs = y[~y_mask]
                X_mis = X[y_mask.data,:] # fall back on numpy bec variablenams
                X_obs = X[~y_mask.data,:]
            

                # A_obs
                propensityModel = LogisticRegression(random_state=0).fit(X,y_mask.data)
                zhat_obs = propensityModel.predict_proba(X_obs)
                A = zhat_obs[:,0] # probability for class 1, ie for missingness
                A_2 = np.square(A)
                A_3 = A_2*A
                A = np.vstack((A,A_2, A_3)).T
                # check if missing values present
                if y_mis.size == 0:
                    logging.info(f'variable {varname} does not have missing values. skip ...')
                    continue

                # regress missing values
                if iverbose:
                    logging.info('regress missing values')
                self.sup.fit(X_obs, y_obs, A_obs)
                y_predict = self.sup.predict(X_mis)

                # update estimates for missing values
                if iverbose:
                    logging.info('update estimates for missing values')
                v = np.where(data.coords['variable'] == varname)[0][0] #TODO more elegant
                data[y_mask,v] = y_predict

            #if self.ispatiotemporal and (self.iter % 400==0): # only every xth iteration
            #    logging.info('update neighborhood features ...')
            #    data = self._update_neighborhood_par(data, varnames)
            #    #data = self._update_neighborhood(data, varnames)
            #    logging.info('update neighborhood features finished')
            #    self.iter_below = 0 # to avoid finishing directly after feature update

            # calculate stopping criterion
            self.delta_n_old = self.delta_n

            # troyanskaya svdimpute: change is below 0.01
            # stekhoven & bühlmann missforest: squared norm difference of imputed values
            # increases for the first time ((mod-obs)**2 / mod**2)
            # self-made convergence heuristic: the denominator is default the standard deviation (i think)
            # for standardised data
            # i.e. algorithm converges if change between steps is smaller than epsilon times std
            self.delta_n = (np.sum(np.abs(data.loc[:,varnames] - data_old.loc[:,varnames])) / 
                            np.mean(np.abs(data.loc[:,varnames]))).item()
            if not iremove:
                maxchange = np.abs(data.loc[:,'skt'] - data_before.loc[:,'skt']).max().item() #REMOVE ONLY FOR IMPUTETHIS
                logging.info(f'new delta: {np.round(self.delta_n, 9)} diff: ' + str(np.round(self.delta_n_old - self.delta_n, 9)) + f' niter: {self.iter} niterbelow: {self.iter_below} max change: ' + str(np.round(maxchange, 8)))
            self.iter += 1 

            if self.iter > self.maxiter:
                self.iter_below += 1

                if self.iter_below > self.n_miniter_below:
                    return data_old

            elif np.abs(self.delta_n_old - self.delta_n) < self.epsilon and self.iter > self.n_miniter:
                self.iter_below += 1

                if self.iter_below > self.n_miniter_below:
                    return data_old
            else:
                self.iter_below = 0


if __name__ == '__main__':

    import xarray as xr
    from sklearn.linear_model import Ridge
    from lose_data import set_random_fraction_to_nan_arr as lose_data
    from sklearn.gaussian_process import GaussianProcessRegressor
    import matplotlib.pyplot as plt
    from anchorRegression import anchorRegression

    # different initialisations for the data
    data = np.arange(100.).reshape(10,10) # linearly increasing: learns perfectly 
    #data = np.random.rand(10,10) # random: learns nothing which is expected
    x,y = np.meshgrid(np.linspace(-1,1,10), np.linspace(-1,1,10))
    d = np.sqrt(x*x+ y*y)
    sigma, mu = 1.,0.
    data = np.exp(-( (d-mu)**2 / (2. * sigma**2))) # gaussian hill
    data = data + np.random.rand(10,10)*0.01 # add some random noise

    data = xr.DataArray(data, coords=[np.arange(10),np.arange(10)], dims=['datapoints','variable']) 

    # lose data
    data_lost = data.copy()
    np.random.seed(0)
    data_lost.values[np.random.rand(*data_lost.shape) < 0.1] = np.nan
    #data_lost = lose_data(data, 0.1)

    data_filled = data_lost.copy()
    data_filled = data_filled.fillna(0)
    #data_filled[np.isnan(data_filled)] = 0
    impute = Imputation(1e-4, anchorRegression, dimension_reduction=None, tradeoff_param = 1, maxiter=100)
    #impute = Imputation(1e-4, GaussianProcessRegressor, dimension_reduction=None, alpha=1, maxiter=100)
    data_filled, _ = impute.impute(data_filled, np.isnan(data_lost), varnames=np.arange(10), iremove=True)
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(30,10))
    ax[0].pcolormesh(data)
    ax[1].pcolormesh(data_filled)
    ax[2].pcolormesh(data - data_filled)
    #plt.pcolormesh(data - data_filled, vmin=-0.01, vmax=0.01)
    plt.show()
