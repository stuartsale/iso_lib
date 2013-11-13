import numpy as np
import math
import scipy.interpolate as si

class iso_grid_tefflogg:
	'''A class which allows isochrones to be queried for many objects at the same time
		- set up to work in {Teff, logg, feh} space
		- Asumes a regular grid in '''

	def __init__(self, filename, metal_col=0, teff_col=1, logg_col=2):
		self.metal_col=metal_col
		self.teff_col=teff_col
		self.logg_col=logg_col


		self.iso_array=np.loadtxt(filename)
		self.metal_dict={}

		self.metal_interp=None

		self.teff_step=0.
		self.teff_min=0.
		self.teff_gridlen=0.

		self.logg_step=0.
		self.logg_min=0.
		self.logg_gridlen=0.

		self.register(self.iso_array)

	def register(self, array):

		print "Registering"

		for i in range(array.shape[0]):
			if array[i,self.metal_col] not in self.metal_dict:
				self.metal_dict[ array[i,self.metal_col] ]=i
#				if self.metal_step==None and i!=0:
#					self.metal_step=i
		print self.metal_dict

		self.metal_interp=si.interp1d(sorted(self.metal_dict.keys()), sorted(self.metal_dict.values()), kind='nearest', bounds_error=False)


		self.teff_min=array[0,self.teff_col]
		for i in range(sorted(self.metal_dict.values())[1]):
			if array[i,self.teff_col]!=self.teff_min:
				self.teff_step=array[i,self.teff_col]-self.teff_min
				self.teff_gridlen=i
				break
		print self.teff_min, self.teff_step, self.teff_gridlen

	
		self.logg_min=array[0,self.logg_col]
		for i in range(self.teff_gridlen):
			if array[i,self.logg_col]!=self.logg_min:
				self.logg_step=array[i,self.logg_col]-self.logg_min
				self.logg_gridlen=i
				break
		print self.logg_min, self.logg_step, self.logg_gridlen
	


	def query(self, feh, teff, logg):

		r=np.zeros(feh.shape)+4.1
		i=np.zeros(feh.shape)+3.9
		ha=np.zeros(feh.shape)+3.81

		return r, i, ha
