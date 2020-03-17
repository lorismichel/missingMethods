"""
convenience functions for era_missing_xr.py

    @author: Verena Bessenbacher
    @date: 08 11 2019
"""

import numpy as np
from scipy.ndimage.filters import generic_filter, median_filter

def temporal_interp(data_one_variable):
    footprint = np.zeros((3,1,1))
    footprint[[0,2],0,0,] = 1
    return generic_filter(data_one_variable, np.nanmean, footprint=footprint, mode='nearest')

def spatial_mean(data_one_variable):
    footprint = np.zeros((1,3,3))
    footprint[0,:,:] = 1
    footprint[0,1,1] = 0
    return generic_filter(data_one_variable, np.nanmean, footprint=footprint, mode='wrap') #TODO: we need wrap for longitude and nearest for latitude

def calc_frac_missing(data, frac_ocean):
    frac_all = np.isnan(data).sum().values / data.size
    frac_land = (frac_all - frac_ocean) / (1. - frac_ocean)
    return frac_land

def mean_impute(data):
    #varmeans = data.mean(dim=('time','longitude','latitude'))
    #varmeans = data.mean(dim=('time','landpoints'))
    # fill each missing spacetimepoint with the temporal average at this spacepoint
    varmeans = data.mean(dim=('time'))
    data_filled = data.fillna(varmeans)
    # this only fills 30% of the missing data because for many spacepoints there are no observation at all
    # second attempt: make a spatial filter to fill with mean of neighboring points
    for v in range(data.shape[0]):
        data_filled[v,:,:] = spatial_mean(data_filled[v,:,:])
    # fill the last remaining holes
    varmeans = data.mean(dim=('time','latitude','longitude'))
    data_filled = data.fillna(varmeans)
    return data_filled

def locf_impute(data):
    data_filled = data.ffill(dim='time')
    # if first value in dim time is nan, nans are retained. fill them with varmeans erase missing values
    #varmeans = data_filled.mean(dim=('time','longitude','latitude'))
    varmeans = data_filled.mean(dim=('time','latitude','longitude'))
    data_filled = data_filled.fillna(varmeans)
    # varmeans method only fills 30% of the missing data because for many spacepoints there are no observation at all
    # third attempt: make a spatial filter to fill with mean of neighboring points
    for v in range(data.shape[0]):
        data_filled[v,:,:] = spatial_mean(data_filled[v,:,:])
    # fill the last remaining holes
    varmeans = data.mean(dim=('time','latitude','longitude'))
    data_filled = data.fillna(varmeans)
    return data_filled

def mse(orig, filled, axis):
    return np.mean((orig - filled)**2, axis=axis)

# long term project: progress bar
# TODO backline does not work
class ProgressBar:

    def __init__(self, ncharacters=100, character='*'):
        self.cumchar = 0
        self.ncharacters = ncharacters
        self.character = character

        print('progress: |')
        self._backline()

    def _backline(self):        
        print('\r', end='') 

    def step(self, fraction):
        nadd = int(fraction * self.ncharacters)
        self.cumchar += nadd
        print(self.character * nadd)

    def finish(self):
        nleft = self.ncharacters - self.cumchar
        self._backline()
        print(self.character * nleft)
        self._backline()
        print('| done')

if __name__ == '__main__':

    pb = ProgressBar()
    pb.step(0.2)
    pb.step(0.4)
    pb.step(0.3)
    pb.finish()
