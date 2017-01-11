from __future__ import print_function, division
import math
import numpy as np
import scipy.interpolate as si

import iso_obj as io


class iso_grid_tefflogg(object):
    ''' Class isogrid
        A class which allows isochrones to be queried.

        This class is designed for isochrones gridded in
        ([M/H], T_eff, log(g)) space.

        Designed so that many queires can be made at
        the same time.
    '''

    def __init__(self, filename, metal_col=0, Mi_col=1, logage_col=2,
                 teff_col=3, logg_col=4, Jac_col=5, bands=None,
                 verbose=False,
                 R_dir="/home/sale/work-oxford/tracks/Phoenix"):
        """ __init__(self, filename, metal_col=0, Mi_col=1,
                     logage_col=2, teff_col=3, logg_col=4,
                     Jac_col=5, bands=None, verbose=False,
                     R_dir="/home/sale/work-oxford/tracks/Phoenix")

            Inititialises an iso_grid_tefflogg object by reading
            in the isochrones from an ascii file.

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
            teff_col : int
                The index of the column in the file that gives
                the effective temperature.
            logg_col : int
                The index of the column in the file that gives
                the surface gravity log(g).
            Jac_col : int
                The index of the column in the file that gives
                the Jacobian for the conversion from
                (initial mass, age) to (T_eff, log(g)).
            bands : list(string)
                The photometric bands which we want to use
            verbose : bool
                If True various information is provided to
                stdout. Useful when debugging.
            R_dir : string
                The location of the files that specify the
                reddening laws.

            Note
            ----
            the band names given in bands must exactly match
            (a subset of) those in the column headings of the
            file and must also match band names used in the
            files that specify the reddening law.
        """

        self.metal_col = metal_col
        self.teff_col = teff_col
        self.logg_col = logg_col

        self.verbose = verbose

        iso_array = np.loadtxt(filename)
        self.metal_dict = {}

        self.metal_interp = None

        self.teff_step = 0.
        self.teff_min = 0.
        self.teff_gridlen = 0.

        self.logg_step = 0.
        self.logg_min = 0.
        self.logg_gridlen = 0.

        self.register(iso_array)

        abs_mags = {}
        with open(filename, 'r') as f:
            first_line = f.readline().split()

        for band in bands:
            abs_mags[band] = iso_array[:, first_line.index(band)-1]

        if metal_col is not None:
            metal_col = first_line.index("[M/H]") - 1
        if Mi_col is not None:
            Mi_col = first_line.index("Mi") - 1
        if logage_col is not None:
            logage_col = first_line.index("logAge") - 1
        if teff_col is not None:
            metal_col = first_line.index("logTe") - 1
        if logg_col is not None:
            metal_col = first_line.index("logg") - 1
        if Jac_col is not None:
            metal_col = first_line.index("Jacobian") - 1

        self.iso_array2 = io.iso_objs(iso_array[:, Mi_col],
                                      iso_array[:, logage_col],
                                      iso_array[:, metal_col],
                                      iso_array[:, teff_col],
                                      iso_array[:, logg_col],
                                      abs_mags, iso_array[:, Jac_col],
                                      bands=bands, R_dir=R_dir)

    def register(self, array):
        """ register(array)

            Establishes the [M/H] values and the steps and
            ranges in T_eff and log(g). These are then
            stored for use when querying the grid.

            Parameters
            ----------
            array : ndarray(float, float)
                The data that will form the isochrone grid, as
                just read in from file.
        """

        if self.verbose:
            print("Registering")

        for i in range(array.shape[0]):

            if array[i, self.metal_col] not in self.metal_dict:
                self.metal_dict[array[i, self.metal_col]] = i

        if self.verbose:
            print("metal dict:", self.metal_dict)

        self.metal_interp = si.interp1d(sorted(self.metal_dict.keys()),
                                        sorted(self.metal_dict.values()),
                                        kind='nearest', bounds_error=False,
                                        fill_value=0)

        self.teff_min = array[0, self.teff_col]
        for i in range(sorted(self.metal_dict.values())[1]):
            if array[i, self.teff_col] != self.teff_min:
                self.teff_step = array[i, self.teff_col]-self.teff_min
                self.teff_gridlen = i
                break
        if self.verbose:
            print("Teff grid:", self.teff_min, self.teff_step,
                  self.teff_gridlen)

        self.logg_min = array[0, self.logg_col]
        for i in range(self.teff_gridlen):
            if array[i, self.logg_col] != self.logg_min:
                self.logg_step = array[i, self.logg_col] - self.logg_min
                self.logg_gridlen = i
                break
        if self.verbose:
            print("logg grid:", self.logg_min, self.logg_step,
                  self.logg_gridlen)

        self.teff_max = np.max(array[:, self.teff_col])
        self.logg_max = np.max(array[:, self.logg_col])

    def query(self, feh, teff, logg):
        """ query( feh, teff, logg)

            query the isochrone grid and obtain the absolute
            magnitudes that correspond to the queried
            ([Fe/H], T_eff, log(g)) .

            Parameters
            ----------
            feh : float, ndarray(float)
                the metallicity of the point(s) to query
            Teff : float, ndarray(float)
                the effective temperature of the point(s)
                to query
            logg : float, ndarray(float)
                the surface gravity of the point(s) to query

            Returns
            -------
            An array of iso_objs objects corresponding to the
            queried values.

            Note
            ----
            The shape of feh, teff and logg must be consistent.
        """

        teff_out_of_bounds = np.logical_or(teff > self.teff_max,
                                           teff < self.teff_min)
        logg_out_of_bounds = np.logical_or(logg > self.logg_max,
                                           logg < self.logg_min)
        out_of_bounds = np.logical_or(teff_out_of_bounds,
                                      logg_out_of_bounds)

        teff[out_of_bounds] = self.teff_max
        logg[out_of_bounds] = self.logg_max

        rows = (self.metal_interp(feh)
                + (np.rint((teff-self.teff_min)/self.teff_step)
                   * self.teff_gridlen)
                + (np.rint((logg-self.logg_min)/self.logg_step)
                   * self.logg_gridlen)).astype(int)

        return self.iso_array2.subset(rows)
