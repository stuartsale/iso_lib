import numpy as np
import scipy.interpolate as si



class R_curves:

    def __init__(self, directory, bands=['r', 'i', 'Ha']):
    
        self.u_splines={}
        self.v_splines={}
        
        for band in bands:
            self.u_splines[band]=[]
            self.v_splines[band]=[]            
        
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
                
      
            R_file=np.genfromtxt("{0:s}/{1:d}.2out".format(directory,int(R*10)), usecols=columns_required_values)[::-1]
               
            for band in bands:
                self.u_splines[band].append(si.UnivariateSpline(np.log10(R_file[:,columns_required_keys.index("Teff")]),
                        R_file[:,columns_required_keys.index("A{}_2".format(band))],k=1) )  
                        
                self.v_splines[band].append(si.UnivariateSpline(np.log10(R_file[:,columns_required_keys.index("Teff")]),
                        R_file[:,columns_required_keys.index("A{}_1".format(band))],k=1) )      
        
