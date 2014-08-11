import numpy as np
import math

from R_curves import R_curves

class iso_objs:
	""" As iso_obj, but intended to hold multiple objects [e.g. an entire library, an isochrone]"""

	def __init__(self, Mi_in, logage_in, feh_in, logT_in, logg_in, r0_in, i0_in, ha0_in, Jac_in=None, 
			ur_in=None, vr_in=None, ui_in=None, vi_in=None, uha_in=None, vha_in=None,
			log_IMF_prob_in=None, log_SFR_prob_in=None, R_set=R_curves("/home/sale/work-oxford/tracks/R_tracks")) :


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
			self.Jac=np.ones(self.Mi.shape)

		if ur_in is None:
			self.ur = np.array([R(self.logT) for R in R_set.ur_splines]).T
		else:
			self.ur = ur_in

		if vr_in is None:
			self.vr = np.array([R(self.logT) for R in R_set.vr_splines]).T
		else:
			self.vr = vr_in

		if ui_in is None:
			self.ui = np.array([R(self.logT) for R in R_set.ui_splines]).T
		else:
			self.ui = ui_in

		if vi_in is None:
			self.vi = np.array([R(self.logT) for R in R_set.vi_splines]).T 
		else:
			self.vi = vi_in

		if uha_in is None:		
			self.uha = np.array([R(self.logT) for R in R_set.uha_splines]).T
		else:
			self.uha = uha_in

		if vha_in is None:
			self.vha = np.array([R(self.logT) for R in R_set.vha_splines]).T
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

	def Ar(self, A0, R=3.1):
		index=np.rint((R-2.1)/0.1)
		return self.ur[index]*A0*A0 + self.vr[index]*A0

	def Ai(self, A0, R=3.1):
		index=np.rint((R-2.1)/0.1)	
		return self.ui[index]*A0*A0 + self.vi[index]*A0

	def Aha(self, A0, R=3.1):
		index=np.rint((R-2.1)/0.1)	
		return self.uha[index]*A0*A0 + self.vha[index]*A0

	def redline(self, r_i1):
		A_int=(-(self.r0-self.i0)+r_i1)/(self.vr-self.vi)
		return (self.ur-self.uha)*pow(A_int,2) + (self.vr-self.vha)*A_int + (self.r0-self.ha0)

	def subset(self, lines):	# Copy certain rows into a new iso_objs instance
		return iso_objs(self.Mi[lines], self.logage[lines], self.feh[lines], self.logT[lines], self.logg[lines], self.r0[lines], self.i0[lines], self.ha0[lines], self.Jac[lines], self.ur[lines], self.vr[lines], self.ui[lines], self.vi[lines], self.uha[lines], self.vha[lines], self.log_IMF_prob[lines], self.log_SFR_prob[lines])


