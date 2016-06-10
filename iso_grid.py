from __future__ import print_function, division
import math
import numpy as np
import scipy.interpolate as si

class iso_grid(object):
    ''' Class isogrid
        A class which allows isochrones to be queried.

        This class is designed for isochrones gridded in 
        ([M/H], initial mass, age) space.

        Designed so that many queires can be made at  
        the same time.
    '''

    def __init__(self, filename, metal_col=0, age_col=1, mass_col=2,
                 regular_age_step=False):
        """ __init__(self, filename, metal_col=0, age_col=1, 
                     mass_col=2, regular_age_step=False)

            Inititialises an iso_grid object by reading in the 
            isochrones from an ascii file.

            Parameters
            ----------
            filename : string
                The file that contains the grid of isochrones
            metal_col : int
                The index of column in the file that contains 
                the metallicity
            age_col : int
                The index of column in the file that contains 
                the age
            mass_col : int
                The index of the column in the file that gives 
                the initial mass.
        """

        # As this class is not complete throw an error
        raise NotImplementedError("This class is only partially "
                                  "written.")

        self.metal_col = metal_col
        self.age_col = age_col
        self.mass_col = mass_col

        self.regular_age_step = regular_age_step

        self.iso_array = np.loadtxt(filename)
        self.metal_dict = {}

        self.metal_interp = None


        if self.regular_age_step:
            self.age_step = 0
            self.age_min = 0
            self.age_gridlen = 0
        else:
            self.age_dict = {}
            self.age_interp = None

        self.mass_dict = {}
        self.mass_interp = None

        self.register(self.iso_array)

    def register(self, array):
        """ register(array)

            Establishes the [M/H] values and the steps and 
            ranges in age and initial mass. These are then 
            stored for use when querying the grid.

            Parameters
            ----------
            array : ndarray(float, float)
                The data that will form the isochrone grid, as
                just read in from file.
        """

        for i in range(array.shape[0]):
            if array[i,self.age_col] not in self.metal_dict:
                self.metal_dict[ array[i,self.metal_col] ] = i


        self.metal_interp = si.interp1d(sorted(self.metal_dict.keys()),
                                     sorted(self.metal_dict.values()),
                                     kind='nearest',
                                     bounds_error=False)

        if self.regular_age_step:
            self.age_min = array[0,self.age_col]
            for i in range(sorted(self.metal_dict.values())[1]):
                if array[i,self.age_col] != self.age_min:
                    self.age_step = array[i,self.age_col]-self.age_min
                    self.age_gridlen = i
                    break
            print(self.age_min, self.age_step, self.age_gridlen)

        else:
            for i in range(sorted(self.metal_dict.values()) ):
                if array[i,self.age_col] not in self.age_dict:
                    self.age_dict[array[i,self.age_col]] = i
            self.age_interp = si.interp1d(sorted(self.age_dict.keys()),
                                       sorted(self.age_dict.values()),
                                       kind='nearest',
                                       bounds_error=False)


        for i in range(self.age_gridlen):
            if array[i,self.mass_col] not in self.mass_dict:
                self.mass_dict[array[i,self.mass_col]] = i


    def query(self, feh, Mi, age):
        """ query( feh, Mi, age)

            query the isochrone grid and obtain the absolute
            magnitudes that correspond to the queried
            ([M/H], initial mass, age) .

            Parameters
            ----------
            feh : float, ndarray(float)
                the metallicity of the point(s) to query
            Mi : float, ndarray(float)
                the initial mass of the point(s) to query
            age : float, ndarray(float)
                the initial mass of the point(s) to query

            Returns
            -------
            r,i,ha : float, ndarray(float)
                absolut magnitudes in the r,i & Ha bands.

            Note
            ----
            The shape of feh, Mi and age must be consistent.
        """

        r = np.zeros(feh.shape)+4.1
        i = np.zeros(feh.shape)+3.9
        ha = np.zeros(feh.shape)+3.81

        return r, i, ha
