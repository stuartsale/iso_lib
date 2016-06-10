import numpy as np
import numpy.ma as ma
import scipy.interpolate as si


class R_curves:
    """ Class R_curves

        Used to describe the extinction laws(s).
        Holds the variables needed for a parametiration of 
        extinction in the form:
        
        A_X = c_X(R)*A_0^2 + d_X(R)*A_0

    """

    def __init__(self, directory, bands):
        """ __init__(directory, bands)
    
            Initialise an R_curves object.

            Parameters
            ----------

            directory : string
                The directory containing the files which
                list the coefficiants that characterise
                the extinction laws
            bands : list(string)
                The photometric bands for which we require
                the reddening laws.
        """
    
    
        self.A1_splines = {}
        self.A2_splines = {}
        
        for band in bands:
            self.A1_splines[band] = []
            self.A2_splines[band] = []            
        
        for R in np.arange(2.1, 5.1, 0.1):
            with open("{0:s}/{1:d}_curves.out".format(directory,
                                              int(R*10)), 'r') as f:
                first_line = f.readline().split()
                
            columns_required_keys = []
            columns_required_values = [] 
                       
            columns_required_keys.append("Teff")
            columns_required_values.append(first_line.index("Teff")-1)
   
             for band in bands:
                columns_required_keys.append("A_{}_1".format(band))
                columns_required_values.append(
                            first_line.index("A_{}_1".format(band))-1)
                
                columns_required_keys.append("A_{}_2".format(band))
                columns_required_values.append(
                            first_line.index("A_{}_2".format(band))-1)
                
      
            R_file = ma.masked_invalid(np.genfromtxt(
                        "{0:s}/{1:d}_curves.out".format(directory,
                                                        int(R*10)),
                        usecols=columns_required_values)[::-1])

            Teff_col = columns_required_keys.index("Teff")   
            Teff_data = np.log10( ma.filled(R_file[:,Teff_col],
                                            fill_value=50000) )

            for band in bands:
                u_col = columns_required_keys.index(
                                                "A_{}_2".format(band))
                v_col = columns_required_keys.index(
                                                "A_{}_1".format(band))  
                

                self.A2_splines[band].append(si.UnivariateSpline(
                    ma.masked_array(Teff_data,
                            mask=R_file[:,u_col].mask).compressed(),
                    R_file[:,u_col].compressed(),k=1) )  
                        
                self.A1_splines[band].append(si.UnivariateSpline(
                    ma.masked_array(Teff_data,
                            mask=R_file[:,v_col].mask).compressed(),
                    R_file[:,v_col].compressed(),k=1) )     
        
