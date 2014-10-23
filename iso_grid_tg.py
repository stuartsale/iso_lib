import numpy as np
import math
import scipy.interpolate as si
import iso_obj as io

class iso_grid_tefflogg:
	'''A class which allows isochrones to be queried for many objects at the same time
		- set up to work in {Teff, logg, feh} space
		- Asumes a regular grid in '''

	def __init__(self, filename, metal_col=0, Mi_col=1, logage_col=2, teff_col=3, logg_col=4, Jac_col=5, bands, verbose=False):
		self.metal_col=metal_col
		self.teff_col=teff_col
		self.logg_col=logg_col
		
		self.verbose=verbose


		iso_array=np.loadtxt(filename)
		self.metal_dict={}

		self.metal_interp=None

		self.teff_step=0.
		self.teff_min=0.
		self.teff_gridlen=0.

		self.logg_step=0.
		self.logg_min=0.
		self.logg_gridlen=0.

		self.register(iso_array)

		abs_mags={}		
		with open(filename, 'r') as f:
			first_line = f.readline().split()
            
		for band in bands:
			abs_mags[band]=iso_array[:,first_line.index(band)-1]
			
		if metal_col is not None:
		    metal_col=first_line.index("[M/H]")-1	
		if Mi_col is not None:
		    Mi_col=first_line.index("Mi")-1
		if logage_col is not None:
		    logage_col=first_line.index("logAge")-1
		if teff_col is not None:
		    metal_col=first_line.index("logTe")-1
		if logg_col is not None:
		    metal_col=first_line.index("logg")-1
		if Jac_col is not None:
		    metal_col=first_line.index("Jacobian")-1		    		    		    		    
		    

		self.iso_array2=io.iso_objs(iso_array[:,Mi_col], iso_array[:,logage_col], iso_array[:,metal_col], iso_array[:,teff_col], iso_array[:,logg_col],
						 abs_mags, iso_array[:,Jac_col])
		

	def register(self, array):

		if self.verbose:
			print "Registering"

		for i in range(array.shape[0]):

			if array[i,self.metal_col] not in self.metal_dict:
				self.metal_dict[ array[i,self.metal_col] ]=i
#				if self.metal_step==None and i!=0:
#					self.metal_step=i
		if self.verbose:
			print "metal dict:", self.metal_dict

		self.metal_interp=si.interp1d(sorted(self.metal_dict.keys()), sorted(self.metal_dict.values()), kind='nearest', bounds_error=False)


		self.teff_min=array[0,self.teff_col]
		for i in range(sorted(self.metal_dict.values())[1]):
			if array[i,self.teff_col]!=self.teff_min:
				self.teff_step=array[i,self.teff_col]-self.teff_min
				self.teff_gridlen=i
				break
		if self.verbose:
			print "Teff grid:", self.teff_min, self.teff_step, self.teff_gridlen

	
		self.logg_min=array[0,self.logg_col]
		for i in range(self.teff_gridlen):
			if array[i,self.logg_col]!=self.logg_min:
				self.logg_step=array[i,self.logg_col]-self.logg_min
				self.logg_gridlen=i
				break
		if self.verbose:
			print "logg grid:", self.logg_min, self.logg_step, self.logg_gridlen
			
		self.teff_max=np.max(array[:,self.teff_col])
		self.logg_max=np.max(array[:,self.logg_col])				
	


	def query(self, feh, teff, logg):
	
		if teff>self.teff_max or teff<self.teff_min or logg>self.logg_max or logg<self.logg_min:
			raise IndexError()

		rows=(self.metal_interp(feh) 
			+ ( np.rint((teff-self.teff_min)/self.teff_step)*self.teff_gridlen ) 
			+ ( np.rint((logg-self.logg_min)/self.logg_step)*self.logg_gridlen )).astype(int)
			
		return self.iso_array2.subset(rows)
		
