import numpy as np
import scipy.interpolate as si



class R_curves:

    def __init__(self, directory):
    
        self.ur_splines=[]
        self.vr_splines=[]
        self.ui_splines=[]
        self.vi_splines=[]
        self.uha_splines=[]
        self.vha_splines=[]
        
        for R in np.arange(2.1, 5.1, 0.1):
            R_file=np.genfromtxt("{0:s}/{1:d}.2out".format(directory,int(R*10)), usecols=(1,13,14,15,16,17,18) )[::-1]
            
            self.ur_splines.append(si.UnivariateSpline(np.log10(R_file[:,0]),R_file[:,2],k=1) )
            self.vr_splines.append(si.UnivariateSpline(np.log10(R_file[:,0]),R_file[:,1],k=1) )
            
            self.ui_splines.append(si.UnivariateSpline(np.log10(R_file[:,0]),R_file[:,4],k=1) )
            self.vi_splines.append(si.UnivariateSpline(np.log10(R_file[:,0]),R_file[:,3],k=1) )
            
            self.uha_splines.append(si.UnivariateSpline(np.log10(R_file[:,0]),R_file[:,6],k=1) )
            self.vha_splines.append(si.UnivariateSpline(np.log10(R_file[:,0]),R_file[:,5],k=1) )
        
        
