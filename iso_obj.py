import numpy as np
import math

class iso_obj:
	""" A class which holds the attributes of a singke object in a set of isochrones"""

	def __init__(self, Mi_in, logage_in, feh_in, logT_in, logg_in, r0_in, i0_in, ha0_in, Jac_in):
		self.Mi=Mi_in
		self.logage=logage_in
		self.feh=feh_in

		self.logT=logT_in
		self.logg=logg_in

		self.r0=r0_in
		self.i0=i0_in
		self.ha0=ha0_in

		self.ur = 0.000537208*(self.r0-self.i0)-0.003828217
		self.vr = −0.030113725*(r0-i0)+0.843291883

		self.ui = 0.000035467*(r0-i0)-0.001929697 
		self.vi = −0.012866951*(r0-i0)+0.602118304 
		
		self.uha = 0.000058201*(r0-i0)-0.000054837
		self.vha = −0.000413472*(r0-i0)+0.762284312 

		self.log_IMF_prob=-2.7*self.Mi
		self.log_SFR_prob=0.		

	def Ar(A0):
		return self.ur*A0*A0 + self.vr*A0

	def Ai(A0):
		return self.ui*A0*A0 + self.vi*A0

	def Aha(A0):
		return self.uha*A0*A0 + self.vha*A0

	def redline(r_i1):
		A_int=(-(self.r0-self.i0)+r_i1)/(self.vr-self.vi)
		return (self.ur-self.uha)*pow(A_int,2) + (self.vr-self.vha)*A_int + (self.r0-self.ha0)

class iso_objs:
	""" As iso_obj, but intended to hold multiple objects [e.g. an entire library, an isochrone]"""

	def __init__(self, Mi_in, logage_in, feh_in, logT_in, logg_in, r0_in, i0_in, ha0_in, Jac_in=None):
		if Mi_in.size!=logage_in or Mi_in.size!=feh_in or Mi_in.size!=logT_in or Mi_in.size!=logg_in or Mi_in.size!=r0_in or Mi_in.size!=i0_in or Mi_in.size!=ha0_in:
			raise TypeError("Mismatched Input Arrays")

		self.Mi=Mi_in
		self.logage=logage_in
		self.feh=feh_in

		self.logT=logT_in
		self.logg=logg_in

		self.r0=r0_in
		self.i0=i0_in
		self.ha0=ha0_in

		if Jac_in:
			self.Jac=Jac_in
		else:
			self.Jac=np.zeros(self.Mi.shape)

		self.ur = 0.000537208*(self.r0-self.i0)-0.003828217
		self.vr = −0.030113725*(r0-i0)+0.843291883

		self.ui = 0.000035467*(r0-i0)-0.001929697 
		self.vi = −0.012866951*(r0-i0)+0.602118304 
		
		self.uha = 0.000058201*(r0-i0)-0.000054837
		self.vha = −0.000413472*(r0-i0)+0.762284312 

		self.log_IMF_prob=-2.7*self.Mi
		self.log_SFR_prob=0.		

	def Ar(A0):
		return self.ur*A0*A0 + self.vr*A0

	def Ai(A0):
		return self.ui*A0*A0 + self.vi*A0

	def Aha(A0):
		return self.uha*A0*A0 + self.vha*A0

	def redline(r_i1):
		A_int=(-(self.r0-self.i0)+r_i1)/(self.vr-self.vi)
		return (self.ur-self.uha)*pow(A_int,2) + (self.vr-self.vha)*A_int + (self.r0-self.ha0)

	def subset(rows):	# Copy certain rows into a new iso_objs instance
