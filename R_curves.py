import numpy as np
import scipy.interpolate as si
import numpy.ma as ma


class R_curves:

    def __init__(self, directory, bands=['r', 'i', 'Ha']):
    
        print "loading R curves"
    
        self.A1_splines={}
        self.A2_splines={}
        
        for band in bands:
            self.A1_splines[band]=[]
            self.A2_splines[band]=[]            
        
        for R in np.arange(2.1, 5.1, 0.1):
            with open("{0:s}/{1:d}.2out".format(directory,int(R*10)), 'r') as f:
                first_line = f.readline().split()
                
            columns_required_keys=[]
            columns_required_values=[] 
                       
            columns_required_keys.append("Teff")
            columns_required_values.append(first_line.index("Teff")-1)
            for band in bands:
                columns_required_keys.append("A{}_1".format(band))
                columns_required_values.append(first_line.index("A{}_1".format(band))-1)
                
                columns_required_keys.append("A{}_2".format(band))
                columns_required_values.append(first_line.index("A{}_2".format(band))-1)
                
      
            R_file=ma.masked_invalid(np.genfromtxt("{0:s}/{1:d}.2out".format(directory,int(R*10)), usecols=columns_required_values)[::-1])

            Teff_col=columns_required_keys.index("Teff")          
            for band in bands:
                u_col=columns_required_keys.index("A{}_2".format(band))
                v_col=columns_required_keys.index("A{}_1".format(band))                

                self.A2_splines[band].append(si.UnivariateSpline(
                    ma.masked_array(np.log10(R_file[:,Teff_col]),mask=R_file[:,u_col].mask).compressed(),
                    R_file[:,u_col].compressed(),k=1) )  
                        
                self.A1_splines[band].append(si.UnivariateSpline(
                    ma.masked_array(np.log10(R_file[:,Teff_col]),mask=R_file[:,v_col].mask).compressed(),
                    R_file[:,v_col].compressed(),k=1) )     
        
