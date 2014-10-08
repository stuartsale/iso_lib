import numpy as np
import math

from R_curves import R_curves

class iso_objs:
	""" As iso_obj, but intended to hold multiple objects [e.g. an entire library, an isochrone]"""

	def __init__(self, Mi_in, logage_in, feh_in, logT_in, logg_in, abs_mag_in,
        	Jac_in=None, AX1_in=None, AX2_in=None,
			log_IMF_prob_in=None, log_SFR_prob_in=None, R_set=R_curves("/home/sale/work-oxford/tracks/Phoenix")) :
			
		bands=abs_mag_in.keys()

#		if Mi_in.size!=logage_in or Mi_in.size!=feh_in or Mi_in.size!=logT_in or Mi_in.size!=logg_in or Mi_in.size!=r0_in or Mi_in.size!=i0_in or Mi_in.size!=ha0_in:
#			raise TypeError("Mismatched Input Arrays")

		self.Mi=Mi_in
		self.logage=logage_in
		self.feh=feh_in

		self.logT=logT_in
		self.logg=logg_in
		
		self.abs_mag=abs_mag_in

		if Jac_in is not None:
			self.Jac=Jac_in
		else:
			self.Jac=np.ones(self.Mi.shape)
			
		if AX1_in is None:
		    self.AX1={}
		    for band in bands:
		        self.AX1[band]=np.array([ R(self.logT) for R in R_set.A1_splines[band] ]).T
		else:
		    self.AX1=AX1_in
		    
		if AX2_in is None:
		    self.AX2={}
		    for band in bands:
		        self.AX2[band]=np.array([ R(self.logT) for R in R_set.A2_splines[band] ]).T
		else:
		    self.AX2=AX2_in		    

		if log_IMF_prob_in is None:
			self.log_IMF_prob=-2.7*self.Mi
		else:
			self.log_IMF_prob=log_IMF_prob_in

		if log_SFR_prob_in is None:
			self.log_SFR_prob=np.zeros(self.Mi.shape)
		else:
			self.log_SFR_prob=log_SFR_prob_in	
			
	def AX(self, band, A0, R=3.1):
		index=np.rint((R-2.1)/0.1)
		return self.AX2[band][index]*A0*A0 + self.AX1[band][index]*A0			
    

	def redline(self, r_i1):
		A_int=(-(self.abs_mag["r"]-self.abs_mag["i"])+r_i1)/(self.Ar1-self.Ai1)
		return (self.Ar2-self.Aha2)*pow(A_int,2) + (self.Ar1-self.Aha2)*A_int + (self.abs_mag["r"]-self.abs_mag["Ha"])

	def subset(self, lines):	# Copy certain rows into a new iso_objs instance
		AX1_cut={}
		AX2_cut={}
		abs_mag_cut={}

		for band in self.abs_mag.keys():
			AX1_cut[band]=self.AX1[band][lines]
			AX2_cut[band]=self.AX2[band][lines]
			abs_mag_cut[band]=self.abs_mag[band][lines]
	    
		return iso_objs(self.Mi[lines], self.logage[lines], self.feh[lines], self.logT[lines], self.logg[lines], abs_mag_cut, self.Jac[lines], 
		        AX1_cut, AX2_cut, self.log_IMF_prob[lines], self.log_SFR_prob[lines])


