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
		self.vr = -0.030113725*(self.r0-self.i0)+0.843291883

		self.ui = 0.000035467*(self.r0-self.i0)-0.001929697 
		self.vi = -0.012866951*(self.r0-self.i0)+0.602118304 
		
		self.uha = 0.000058201*(self.r0-self.i0)-0.000054837
		self.vha = -0.000413472*(self.r0-self.i0)+0.762284312 

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

	def __init__(self, Mi_in, logage_in, feh_in, logT_in, logg_in, r0_in, i0_in, ha0_in, Jac_in=None, 
			ur_in=None, vr_in=None, ui_in=None, vi_in=None, uha_in=None, vha_in=None,
			log_IMF_prob_in=None, log_SFR_prob_in=None) :


#		if Mi_in.size!=logage_in or Mi_in.size!=feh_in or Mi_in.size!=logT_in or Mi_in.size!=logg_in or Mi_in.size!=r0_in or Mi_in.size!=i0_in or Mi_in.size!=ha0_in:
#			raise TypeError("Mismatched Input Arrays")

		self.Mi=Mi_in
		self.logage=logage_in
		self.feh=feh_in

		self.logT=logT_in
		self.logg=logg_in

		self.r0=r0_in
		self.i0=i0_in
		self.ha0=ha0_in

		if Jac_in is not None:
			self.Jac=Jac_in
		else:
			self.Jac=np.zeros(self.Mi.shape)

		if ur_in is None:
			self.ur = 0.000537208*(self.r0-self.i0)-0.003828217
		else:
			self.ur = ur_in

		if vr_in is None:
			self.vr = -0.030113725*(self.r0-self.i0)+0.843291883
		else:
			self.vr = vr_in

		if vr_in is None:
			self.ui = 0.000035467*(self.r0-self.i0)-0.001929697 
		else:
			self.ui = ui_in

		if vr_in is None:
			self.vi = -0.012866951*(self.r0-self.i0)+0.602118304 
		else:
			self.vi = vi_in

		if vr_in is None:		
			self.uha = 0.000058201*(self.r0-self.i0)-0.000054837
		else:
			self.uha = uha_in

		if vr_in is None:
			self.vha = -0.000413472*(self.r0-self.i0)+0.762284312
		else:
			self.vha = vha_in

		if log_IMF_prob_in is None:
			self.log_IMF_prob=-2.7*self.Mi
		else:
			self.log_IMF_prob=log_IMF_prob_in

		if log_SFR_prob_in is None:
			self.log_SFR_prob=np.zeros(self.Mi.shape)
		else:
			self.log_SFR_prob=log_SFR_prob_in	

	def Ar(self, A0):
		return self.ur*A0*A0 + self.vr*A0

	def Ai(self, A0):
		return self.ui*A0*A0 + self.vi*A0

	def Aha(self, A0):
		return self.uha*A0*A0 + self.vha*A0

	def redline(self, r_i1):
		A_int=(-(self.r0-self.i0)+r_i1)/(self.vr-self.vi)
		return (self.ur-self.uha)*pow(A_int,2) + (self.vr-self.vha)*A_int + (self.r0-self.ha0)

	def subset(self, lines):	# Copy certain rows into a new iso_objs instance
		return iso_objs(self.Mi[lines], self.logage[lines], self.feh[lines], self.logT[lines], self.logg[lines], self.r0[lines], self.i0[lines], self.ha0[lines], self.Jac[lines], self.ur[lines], self.vr[lines], self.ui[lines], self.vi[lines], self.uha[lines], self.vha[lines], self.log_IMF_prob[lines], self.log_SFR_prob[lines])


