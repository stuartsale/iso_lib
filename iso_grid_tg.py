import numpy as np
import math
import scipy.interpolate as si
import iso_obj as io

class iso_grid_tefflogg:
	'''A class which allows isochrones to be queried for many objects at the same time
		- set up to work in {Teff, logg, feh} space
		- Asumes a regular grid in '''

	def __init__(self, filename, metal_col=0, teff_col=1, logg_col=2, verbose=False):
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

		self.iso_array2=io.iso_objs(iso_array[:,1], iso_array[:,2], iso_array[:,0], iso_array[:,3], iso_array[:,4],
						 iso_array[:,9], iso_array[:,10], iso_array[:,11], iso_array[:,5])
		

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
	


	def query(self, feh, teff, logg):

		rows=(self.metal_interp(feh) 
			+ ( np.rint((teff-self.teff_min)/self.teff_step)*self.teff_gridlen ) 
			+ ( np.rint((logg-self.logg_min)/self.logg_step)*self.logg_gridlen )).astype(int)

#		r=self.iso_array[rows, 9]
#		i=self.iso_array[rows, 10]
#		ha=self.iso_array[rows, 11]

#		return r, i, ha
		return self.iso_array2.subset(rows)
		
