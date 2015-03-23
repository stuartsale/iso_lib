import numpy as np
import math
import scipy.interpolate as si

class iso_grid:
    '''A class which allows isochrones to be queried for many objects at the same time'''

    def __init__(self, filename, metal_col=0, age_col=1, mass_col=2, regular_age_step=False):
        self.metal_col=metal_col
        self.age_col=age_col
        self.mass_col=mass_col

        self.regular_age_step=regular_age_step

        self.iso_array=np.loadtxt(filename)
        self.metal_dict={}

        self.metal_interp=None


        if self.regular_age_step:
            self.age_step=0
            self.age_min=0
            self.age_gridlen=0
        else:
            self.age_dict={}
            self.age_interp=None

        self.mass_dict={}
        self.mass_interp=None

        self.register(self.iso_array)

    def register(self, array):

        for i in range(array.shape[0]):
            if array[i,self.age_col] not in self.metal_dict:
                self.metal_dict[ array[i,self.metal_col] ]=i


        self.metal_interp=si.interp1d(sorted(self.metal_dict.keys()), sorted(self.metal_dict.values()), kind='nearest', bounds_error=False)

        if self.regular_age_step:
            self.age_min=array[0,self.age_col]
            for i in range(sorted(self.metal_dict.values())[1]):
                if array[i,self.age_col]!=self.age_min:
                    self.age_step=array[i,self.age_col]-self.age_min
                    self.age_gridlen=i
                    break
            print self.age_min, self.age_step, self.age_gridlen

        else:
            for i in range(sorted(self.metal_dict.values()) ):
                if array[i,self.age_col] not in self.age_dict:
                    self.age_dict[ array[i,self.age_col] ]=i
            self.age_interp=si.interp1d(sorted(self.age_dict.keys()), sorted(self.age_dict.values()), kind='nearest', bounds_error=False)


        for i in range(self.age_gridlen):
            if array[i,self.mass_col] not in self.mass_dict:
                self.mass_dict[ array[i,self.mass_col] ]=i


    def query(self, feh, Mi, age):

        r=np.zeros(feh.shape)+4.1
        i=np.zeros(feh.shape)+3.9
        ha=np.zeros(feh.shape)+3.81

        return r, i, ha
