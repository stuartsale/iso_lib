import numpy as np
import math

from R_curves import R_curves

class iso_objs:
    """ A class designed to hold an array of points from
        (a set of) isochrone(s).

        Each instance could hold an entire library, a 
        single isochrone or a set of random objects.
    """

    def __init__(self, Mi_in, logage_in, feh_in, logT_in, logg_in, abs_mag_in,
            Jac_in=None, AX1_in=None, AX2_in=None,
            log_IMF_prob_in=None, log_SFR_prob_in=None, R_set=None, bands=None,
            R_dir="/home/sale/work-oxford/tracks/Phoenix") :#R_curves("/home/sale/work-oxford/tracks/Phoenix")
        """ __init__(self, Mi_in, logage_in, feh_in, logT_in, 
                     logg_in, abs_mag_in, Jac_in=None, 
                     AX1_in=None, AX2_in=None, 
                     log_IMF_prob_in=None, log_SFR_prob_in=None, 
                     R_set=None, bands=None,
                     R_dir="/home/sale/work-oxford/tracks/Phoenix")

            Initialise an iso_objs object. 
            Not intended for direct use by a user, should 
            be called by other functions.
            
            Parameters
            ----------
            Mi_in : ndarray(float)
                The initial masses of the objects
            log_in : ndarray(float)
                The logarithm of the ages of the objects
            feh_in : ndarray(float)
                The metallicities of the objects
            logT_in : ndarray(float)
                The log(T_eff) of the objects
            logg_in : ndarray(float)
                The log(g) of the objects
            abs_mag_in : dict(string:ndarray(float))
                The absolute magnitudes of the objects
            Jac_in : ndarray(float), optional
                The Jacobian of the (initial mass, age) 
                to (T_eff, log(g)) conversion for the objects.
                If not given they are all assumed to be 1.
            AX1_mag_in : dict(string:float), optional
                The linear coefficiants in the parametrisation
                of extinction the objects as 
                A_X = c*A_0^2 + d*A_0 .
                If not given they are found from the given
                set of R curves.
            AX1_mag_in : dict(string:float), optional
                The quadratic coefficiants in the 
                parametrisation of extinction the objects as 
                A_X = c*A_0^2 + d*A_0.
                If not given they are found from the given
                set of R curves.
            log_IMF_prob_in : ndarray(float), optional
                The logarithms of the IMF probs of the objects
                given their initial masses.
                If not given they are calculated assuming
                p(M_i) \propto M_i^{-2.7} .
            log_SFR_prob_in : ndarray(float), optional
                The logarithms of the SFR probs of the objects
                given their initial masses.
                If not given they are calculated assuming
                a uniform star formation rate.
            R_set : R_curves, optional
                A set of reddening laws
            R_dir : string
                The location of the files that specify the 
                reddening laws.    
        """
            
        if  AX1_in is None and AX2_in is None:
            if bands is not None:
                R_set=R_curves(R_dir, bands)
            else:            
                raise ValueError("Either the bands needed must be given or extinction coefficiants provided")

        self.Mi=Mi_in
        self.logage=logage_in
        self.feh=feh_in

        self.logT=logT_in
        self.logg=logg_in
        
        self.abs_mag=abs_mag_in

        if Jac_in is not None:
            self.Jac=Jac_in
        else:
            self.Jac=np.ones(self.Mi.shape)
            
        if AX1_in is None:
            self.AX1={}
#            print "splines: ", R_set.A1_splines.keys()
            for band in bands:
                self.AX1[band]=np.array([ R(self.logT) for R in R_set.A1_splines[band] ]).T
        else:
            self.AX1=AX1_in
            
        if AX2_in is None:
            self.AX2={}
            for band in bands:
                self.AX2[band]=np.array([ R(self.logT) for R in R_set.A2_splines[band] ]).T
        else:
            self.AX2=AX2_in            

        if log_IMF_prob_in is None:
            self.log_IMF_prob=-2.7*self.Mi
        else:
            self.log_IMF_prob=log_IMF_prob_in

        if log_SFR_prob_in is None:
            self.log_SFR_prob=np.zeros(self.Mi.shape)
        else:
            self.log_SFR_prob=log_SFR_prob_in    
            
    def AX(self, band, A0, R=3.1):
        """ AX(band, A0, R=3.1)
        
            Get the redding in some band, for a set of objects
            assuming A_X = c*A_0^2 + d*A_0.

            Parameters
            ----------
            band : string, ndarray(string)
                The photometric band(s) for which we 
                require extinctions. 
            A0 : float, ndarray(float)
                Monochromatic extinction.
                Either a single float in which case the 
                band extinction for all objects subject to
                this monochromatic extinction.
                Or an ndarray of floats of the same shape as 
                the contents of iso_objs, in which case this
                is the monochromatic extinctions for each 
                object.
            R : float, ndarray(float), optional
                The parameter R_5495 that specifies the 
                shape of the reddening law to be used.
                As with A0, can be either the same value for 
                all isochrone points or an array of values so 
                that each point gets its own value.

            Returns
            -------
            A_X : ndarray(float)
                The band extinction
        """

        index=np.rint((R-2.1)/0.1).astype('int')
        return self.AX2[band][np.arange(index.size),index]*A0*A0 + self.AX1[band][np.arange(index.size),index]*A0            
    

    def redline(self, r_i1):
        """ redline(r_i1)
            
            Provides a reddening line in IPHAS colour space.

            Pararameters
            ------------
            r_i1 : float, ndarray(float)
                the (r-i) colour(s) to be queried.
            
            Returns
            -------
            float, ndarray(float)
                The corresponding (r-Ha) colour(s)
        """
        A_int=(-(self.abs_mag["r"]-self.abs_mag["i"])+r_i1)/(self.Ar1-self.Ai1)
        return (self.Ar2-self.Aha2)*pow(A_int,2) + (self.Ar1-self.Aha2)*A_int + (self.abs_mag["r"]-self.abs_mag["Ha"])

    def subset(self, lines):
        """ subset(lines)
        
            Copy certain rows into a new iso_objs instance

            Parameters
            ----------
            lines : ndarray(int)
                The indices of the objects to be copied

            Returns
            -------
            iso_objs
                An iso_objs object containing the copied 
                isochrone points.
        """

        AX1_cut={}
        AX2_cut={}
        abs_mag_cut={}

        for band in self.abs_mag.keys():
            AX1_cut[band]=self.AX1[band][lines]
            AX2_cut[band]=self.AX2[band][lines]
            abs_mag_cut[band]=self.abs_mag[band][lines]
        
        return iso_objs(self.Mi[lines], self.logage[lines], self.feh[lines], self.logT[lines], self.logg[lines], abs_mag_cut, self.Jac[lines], 
                AX1_cut, AX2_cut, self.log_IMF_prob[lines], self.log_SFR_prob[lines])


