# This contains all the function that describes the Bumped period spacing
# and how they can be used to derived mode amplitudes from inertia ratio
# as they have been developed and tested. This arises from the following publications:
# https://arxiv.org/pdf/1509.06193.pdf (Inertia and ksi relation)
# https://www.aanda.org/articles/aa/pdf/2015/08/aa26449-15.pdf (Eq. 17 for the rotation - splitting relation)
# https://arxiv.org/pdf/1401.3096.pdf (Fig 13 and 14 for the evolution - rotation relation in SG and RGB) 
# https://arxiv.org/pdf/1505.06087.pdf (Eq. 3 for determining log(g) using Teff and numax, used for getting the evolution stage in conjonction with Fig. 13 from above)
# https://iopscience.iop.org/article/10.1088/2041-8205/781/2/L29/pdf (Inertia and Height relation)
# https://arxiv.org/pdf/1707.05989.pdf (Fig.5, distribution of surface rotation for 361 RGB stars)
import numpy
import matplotlib.pyplot as plt
from scipy import interpolate
import solver_mm


# the ksi function as defined in equation 14 of Mosser+2017 (https://arxiv.org/pdf/1509.06193.pdf)
# Inputs: 
# 	nu: The freqency axis (microHz)
# 	nu_p : Frequenc of a p-mode (microHz)
#	nu_g : Frequency of a g-mode (microHz)
#	Dnu_p: Large separation. DOES NOT NEED TO BE A CONSTANT (fonction of np or nu to account for glitches)
#   DPl: Period Spacing (seconds)
#	q : Coupling term (no unit)
def ksi_fct1(nu, nu_p, nu_g, Dnu_p, DPl, q):
	
	cos_upterm=numpy.pi * 1e6 * (1./nu - 1./nu_g)/DPl
	cos_downterm=numpy.pi * (nu - nu_p) /Dnu_p
	front_term= 1e-6 * nu**2 * DPl / (q * Dnu_p) # relation accounting for units in Hz and in seconds

	ksi=1./(1. + front_term * numpy.cos(cos_upterm)**2/numpy.cos(cos_downterm)**2)
	
	return ksi

# Variant of ksi_fct that deals with arrays for nu_p, nu_g, Dnu_p, DPl
# This requires that nu_p and Dnu_p have the same dimension
# Also nu_g and DPl must have the same dimension
# Additional parameter:
#  - norm-method: When set to 'fast', normalise by the max of the ksi_pg calculated at the
#				   Frequencies nu given by the user
#				   When set to 'exact', normalise by the max of a heavily interpolated function
#				   of ksi_pg. Allows a much a higher precision, but will be slower
#				   This could be usefull as in case of low ng, the norm is badly estimated in
#				   'fast' mode. Then we need to use a more continuous function to evaluate the norm
def ksi_fct2(nu, nu_p, nu_g, Dnu_p, DPl, q, norm_method='fast'):
	Lp=len(nu_p)
	Lg=len(nu_g)

	ksi_pg=0.
	for np in range(Lp):
		for ng in range(Lg):
			ksi_tmp=ksi_fct1(nu, nu_p[np], nu_g[ng], Dnu_p[np], DPl[ng], q)
			ksi_pg=ksi_pg + ksi_tmp

	if norm_method == 'fast':
		norm_coef=max(ksi_pg)
	else: # We build a very resolved 'continuous function of the frequency to calculate the norm'
		fmin=min([min(nu_p), min(nu_g)])
		fmax=max([max(nu_p), max(nu_g)])
		resol=0.1
		Ndata=int((fmax-fmin)/resol)
		nu4norm=numpy.linspace(fmin, fmax, Ndata)
		ksi4norm=0.
		for np in range(Lp):
			for ng in range(Lg):
				ksi_tmp=ksi_fct1(nu4norm, nu_p[np], nu_g[ng], Dnu_p[np], DPl[ng], q)
				ksi4norm=ksi4norm + ksi_tmp
		norm_coef=max(ksi4norm)

	ksi_pg=ksi_pg/norm_coef	
	return ksi_pg

# Inertia as theoretically defined in function of ksi
def inertia_ratio_ksi(ksi):
	IlI0=1./(1. - ksi)
	return IlI0

# Inertia as theoretically defined in function of the amplitude and the widths of the modes
# This could be used to approximate the inertia of the modes from observations of amplitudes and widths
# see https://iopscience.iop.org/article/10.1088/2041-8205/781/2/L29/pdf for further info 
# Inputs:
#	amp0: Amplitudes of the l=0 modes (ppm)
#   ampl: Amplitudes of the l mixed modes (ppm)
#   visibility_l_square: Mode visibility defined as Vl^2 = H1 / H0. Ex. V(l=1)^2 = 1.5
#   width0 : Widths of the l=0 modes
#	widthl : Widths of the l mixed modes
# Outputs:
#	IlI0: Inertia ratio
def inertia_ratio_obs_amp(amp0, ampl, visibility_l_square, width0, widthl):
	r1=amp0/ampl
	r2=numpy.sqrt(width0/widthl)
	IlI0=numpy.sqrt(visibility_l)*r1*r2
	return IlI0

# Inertia as theoretically defined in function of the height and the widths of the modes
# This could be used to approximate the inertia of the modes from observations of amplitudes and widths
# see https://iopscience.iop.org/article/10.1088/2041-8205/781/2/L29/pdf for further info 
# Inputs:
#	h0: Heights of the l=0 modes (ppm)
#   hl: Heights of the l mixed modes (ppm)
#   visibility_l_square: Mode visibility defined as Vl^2 = H1 / H0. Ex. V(l=1)^2 = 1.5
#   width0 : Widths of the l=0 modes
#	widthl : Widths of the l mixed modes
# Outputs:
#	IlI0: Inertia ratio
def inertia_ratio_obs_height(h0, hl, visibility_l_square, width0, widthl):
	r1=numpy.sqrt(h0/h1)
	r2=width0/width1
	IlI0=visibility_l_square*r1*r2
	return IlI0

# Original function is comming from the IDL function 'generate_cfg_asymptotic_Hgauss_Wscaled_act_asym'
# available in /Volumes/MCMC_RES/Repository/Spectra-Simulator-IDL/v1.4/
# This width profile obviously applies only to l=0 modes OR to any l IF the star is in the main sequence.
# Inputs:
#	- nu_star: Frequencies for the star
#	- Dnu_star: The large separation
#	- n_at_numax_star: The 
# Outputs:
#	- w: A rescaled Solar width profile, normalised to 1 at numax. To be used to make a synthetic star mimmicking the Solar profile
#	- h: A rescaled Solar height profile, normalised to 1 at numax. To be used to make a synthetic star mimmicking the Solar profile
def width_height_MS_sun_rescaled(nu_star, Dnu_star, numax_star):
	# ***** CONSTANTS *****
	Dnu_sun=135.1
	epsilon_sun=0.5
	#numax_sun=3150. # THIS IS THE NUMAX IN AMPLITUDE
	numax_sun=2900. # THIS IS THE NUMAX IN HEIGHT

	# -------- Relation giving the solar width profile from Appourchaux et al. 2014 ------
	alfa=4.97
	Gamma_alfa=4.65
	Wdip=4646. # depth of the dip
	nu_dip=3083.
	Ddip=4.66

	# Frequencies of the l=0 of the Sun + [1650, 1750, 1850] + [4180, 4280, 4380, 4480] to avoid extrapolation and favor interpolation
	nu_sun=[1957.4748, 2093.5983, 2228.8442, 2362.8797, 2496.3328, 2629.8436, 2764.3597, 2899.2249, 3033.9623, 3168.9156, 3303.8225, 3439.3876, 3575.2118, 3711.6045, 3848.5361, 3984.6612]
	nu_sun=numpy.asarray(nu_sun)
	
	height_sun=[0.55633623, 0.71080326, 0.84916942, 1.0309479, 1.3676815, 2.0930273, 2.8720608, 3.9032770, 3.7507970, 2.8629352, 1.8167902, 0.92533429, 0.42467669, 0.17490098, 0.079882521, 0.038872344]
	height_sun=numpy.asarray(height_sun)
	int_fct_hsun = interpolate.interp1d(nu_sun, height_sun)
	height_sun_at_numax=int_fct_hsun(numax_sun)

	lnGamma0= alfa*numpy.log(nu_sun/numax_sun) +  numpy.log(Gamma_alfa)
	lnLorentz= -numpy.log( Ddip )/ (1. + (2. *numpy.log(nu_sun/nu_dip)/numpy.log(Wdip/numax_sun))**2)
	Gamma_sun=numpy.exp(lnGamma0 + lnLorentz)
	int_fct_wsun = interpolate.interp1d(nu_sun, Gamma_sun)
	Gamma_sun_at_numax=int_fct_wsun(numax_sun)

	n_at_numax_sun=numax_sun/Dnu_sun - epsilon_sun 
	en_list_sun=nu_sun/Dnu_sun - epsilon_sun # This list will be monotonic from 14 until 29 (with step ~1)
	# ------------------------------------------------------------------------------------

	# Rescaling using the base frequencies given above for the Sun
	epsilon_star=nu_star/Dnu_star % 1 # NEED SOME CHECK REGARDING EPSILON CALCULATION
	n_at_numax_star=numax_star/Dnu_star - epsilon_star

	int_fct_w = interpolate.interp1d(en_list_sun - n_at_numax_sun, Gamma_sun/Gamma_sun_at_numax)
	w=int_fct_w(nu_star/Dnu_star - n_at_numax_star)

	int_fct_h = interpolate.interp1d(en_list_sun - n_at_numax_sun, height_sun/height_sun_at_numax)
	h=int_fct_h(nu_star/Dnu_star - n_at_numax_star)

	return w,h

# Using Inertias, heights and widths of l=0 modes, 
# we can in principe know what is the width of the l mixed modes
# BEWARE HERE: I assume hl/h0(nu) = constant (1 by default). It is not certain
# that this is really supported by observations and theory. 
# Also, this relation does not account for bolometric correction (hl/h0=1.5)
# Inputs:
#	el (integer): degree of the modes
#	nu_m  (array): frequencies for the mixed modes 
#	nu_p (array): frequencies for the p modes
#   nu_g (array): frequencies for the g modes
#	Dnu_p (array): Large separation in function of the frequency of the l=el p modes
#	DPl (array): Period spacing in function of the frequency of the l=el g modes
#	q (double): Coupling coefficient
#	nu_p0 (array):  Frequencies of the l=0 modes. Used for the interplation of the maximum width of el modes
#	width0 (array): Width of l=0 modes. Must be of same size as nu_p0. This array is interpolated
#		at the frequencies of l=el mixed modes. Extrapolation is forbidden (expect a crash)
def gamma_l_fct1(nu_m, nu_p, nu_g, Dnu_p, DP1, q, nu_p0, width0, el, hl_h0_ratio=1):
	#Lp=len(nu_p)
	#Lg=len(nu_g)
	ksi_pg=ksi_fct2(nu_m, nu_p, nu_g, Dnu_p, DP1, q)

	if len(nu_p0) != len(width0):
		print('Inconsistency between the size of the Width and l=0 frequency array')
		print('Cannot pursue. The program will exit now')
		exit()
	else:
		# Perform the interpolation
		int_fct = interpolate.interp1d(nu_p0, width0)
		width0_at_l=int_fct(nu_m)
		width_l=width0_at_l * (1. - ksi_pg)/ numpy.sqrt(hl_h0_ratio)

	return width_l

def gamma_l_fct2(ksi_pg, nu_m, nu_p0, width0, hl_h0_ratio, el):

	if len(nu_p0) != len(width0):
		print('Inconsistency between the size of the Width and l=0 frequency array')
		print('Cannot pursue. The program will exit now')
		exit()
	else:
		# Perform the interpolation
		int_fct = interpolate.interp1d(nu_p0, width0)
		width0_at_l=int_fct(nu_m)
		width_l=width0_at_l * (1. - ksi_pg)/ numpy.sqrt(hl_h0_ratio)
		#print('ksi_pg     ,   width_l    ,   hl_h0_ratio')
		#for i in range(len(ksi_pg)):
		#	print(ksi_pg[i],    width_l[i],   hl_h0_ratio[i])
	return width_l

def h_l_rgb(ksi_pg):
	hl_h0=numpy.sqrt(1. - ksi_pg)
	hl_h0[hl_h0 == 0] = 1e-10
	#for i in range(len(ksi_pg)):
	#	print(ksi_pg[i], hl_h0[i])
	return hl_h0

# Splitting of modes assuming a two-zone (or two-zone averaged) rotation profile
# rot_envelope and rot_core must be in Hz units (or micro, nano, etc...)
# ksi_pg: The ksi function that describes the degree of mixture between p and g modes in function of the more frequency
# rot_envelope: average rotation in the envelope. Must be a scalar
# rot_core: average rotation in the core. Must be a scalar
# Returns:
#	dnu_rot: A vector of same size as ksi_pg
def dnu_rot_2zones(ksi_pg, rot_envelope, rot_core):
	#dnu_rot=(ksi_pg*(rot_core/2 - rot_envelope) + rot_envelope)/(2.*numpy.pi)
	dnu_rot=ksi_pg*(rot_core/2 - rot_envelope) + rot_envelope
	return dnu_rot

# Function that determine the rotation in the envelope, here approximated to be the surface rotation.
# Inspired by the surface rotation from Ceillier et al. 2017 (https://arxiv.org/pdf/1707.05989.pdf), Fig. 5
# They have a skewed distribution going for ~30 days to ~160 days with a peak around 60 days. 
# For simplicity, I just insert a truncated gaussian distribution with rotation between 30  and 90 and a median of 60.
# The truncation happens at sigma
# Returns: 
#	rot_s: rotation frequency in microHz
def rot_envelope():
	sigma=3.
	med=60. # in days
	var=30./sigma # in days
	period_s=numpy.random.normal(loc=med, scale=var)
	if period_s < med - var*sigma:
		period_s=med-var*sigma
	if period_s > med + var*sigma:
		period_s=med+var*sigma

	rot_s=1e6/(86400.*period_s)

	return rot_s

# Use various relations from different papers for getting core rotation frequency
# in fonction of the evolutionary stage of the RGB, and using log(g) as an age proxy
# This can then be used by dnu_rot_2zones() to determine the rotational splittings
# rot_envelope: The rotation of the envelope must be provided (in Hz units)
# numax: frequency at maximum height of the modes of the star
# Teff: Effective temperature of the star. For a RGB/SG, typically between 5500K and 4000K. 
# sigma_logg: Uncertainty to be used on loggg if randomize==True
# randomize_logg: If True, add some random gaussian noise to logg, as per defined in sigma_logg (absolute error)
# randomize_core_rot: If True, add some random gaussian noise to the rotation rate of the core, as per derived by observations from Deheuvels 2014 (see below).
# Returns: 
# 	 (1) rot_envelope: average rotation in the envelope
#	 (2) rot_core: average rotation in the core 
#	 (3) rot_envelope_err: 1sigma error (or dispersion) on rot_envelope 
#	 (4) rot_core_err: 1sigma error (or dispersion) on rot_core
def rot_2zones_v1(rot_envelope, numax, Teff, sigma_logg=0.1, randomize_logg=True, randomize_core_rot=True):
	# Get log(g) from T and numax using https://arxiv.org/pdf/1505.06087.pdf
	numax_sun=2900. # THIS IS THE NUMAX IN HEIGHT
	Teff_sun=5777.
	gsun=28020 # cgs
	g=(numax/numax_sun)*numpy.sqrt(Teff/Teff_sun) * gsun
	logg=numpy.log10(g)
	#print('log(g_star):', logg)
	if randomize_logg == True:
		logg=numpy.random.normal(loc=logg, scale=sigma_logg)
	#print('log(g_star) (randomized):', logg)
	# If an SG/RGB Use the best fit for the core/envelope contrast established using the relation from Fig. 13 of https://arxiv.org/pdf/1401.3096.pdf
	# The used table for the fit is the following:
	# star    / rot ratio (OLA)  / log(g) (graphic estimate)
	# A          2.5 +/- 0.7       3.823 +/- 0.0375
	# B          3.8 +/- 1.1       3.77  +/- 0.02
	# C          6.9 +/ 1          3.76 +/- 0.04 
	# D          15.1 +/- 4        3.71 +/- 0.03
	# E          10.7 +/- 1.5      3.68 +/- 0.02
	# F          21   +/- 8.2      3.6  +/- 0.02
	# using polyfit (only y-uncertainty accounted for), we get:
	# ratio(logg) = a  + b log(g)  + c log(g)^2 
	#a = 1797.04 (+/- 2583.32)   /   b= -891.294 +/- 1374.98    /  c= 110.351  +/- 182.934
	# Uncertainty on each data point from the fit :  2.67476      2.62932      2.79434      3.39917      5.02476      17.3518
	# Best curve on each data point               :  2.43628      5.26652      5.87004      9.21873      11.4928      18.5280
	# Ratio of the two                            :  1.09789     0.499253     0.476035     0.368724     0.437210     0.936516
	# Mean relative quadratic uncertainty                   : 0.693  (69.3%)
	# Mean quadratic uncertainty without the two extremes    : 0.448  (44.8%)

	a=1797.04
	b=-891.294
	c=110.351

	logg_min=3.55
	logg_max=3.8825
	if logg >= logg_min and logg < logg_max : # Case excluding old rgb and main sequence stars. 3.8825 allows slightly faster envelope vs core, but it is tolerable
		core2envelope_star=a + b*logg + c*logg**2
		core2envelope_err_star= core2envelope_star * 0.448 # use the average relative uncertainty calculated above without the extreme data points
	# If an old RGB constant core/envelope contrast (ARBITRARY AS THE ROTATION IS NORMALY STILL SLOWING DOWN IN THE ENVELOPE)
	if logg < logg_min:
		core2envelope_star=a + b*logg_min + c*logg_min**2
		core2envelope_err_star= core2envelope_star * 0.448 # use the average relative uncertainty calculated above without the extreme data points
	# If a MS star, assume uniform rotation (core/envelope == 1)
	if logg >= logg_max:
		core2envelope_star=1.
		core2envelope_err_star=0.1 # arbitrary error for MS case

	#print('core2envelope_star:', core2envelope_star, '  +/-', core2envelope_err_star)

	# Determine the core rotation
	rot_core=core2envelope_star * rot_envelope

	rot_envelope_true=rot_envelope
	rot_core_true=rot_core

	#print('True core rotation: ', rot_core_true)
	#print('True env rotation: ', rot_envelope_true)
	if randomize_core_rot == True:
		# We will assume that the error contribution from the core and from the envelope is the same in relative value:
		# err(rot_core)/rot_core = err(rot_envelope)/rot_envelope ===> core2envelope_err_star = X * sqrt(1 + core2envelope_star**2) 
		# with X = err(rot_core)/rot_core == err(rot_envelope)/rot_envelope... X is what we need here.
		rot_core_err=core2envelope_err_star * rot_core / numpy.sqrt(1 + core2envelope_star**2)
		rot_envelope_err=core2envelope_err_star * rot_envelope / numpy.sqrt(1 + core2envelope_star**2)

		rot_envelope=numpy.random.normal(loc=rot_envelope, scale=rot_envelope_err)
		rot_core=numpy.random.normal(loc=rot_core, scale=rot_core_err)
	else:
		rot_core_err=0
		rot_envelope_err=0

	#print(' numpy.sqrt(1 + core2envelope_star**2) = ', numpy.sqrt(1 + core2envelope_star**2))
	#print('rot_envelope_err=', rot_envelope_err)
	#print('rot_core_err=', rot_core_err)
	#print('rot_envelope (randomized): ', rot_envelope)
	#print('rot_core (randomized): ', rot_core)
	#exit()
	return rot_envelope, rot_core, rot_envelope_true, rot_core_true

#Assumptions: nu_max is derived from the requested Dnu_star parameter using the relation from Stello+2009. 
#	Dnu ~ 0.263*numax^0.77 (no uncertainty implemented here)
def numax_from_stello2009(Dnu_star):
	# Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	beta0=0.263 # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	beta1=0.77 # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	nu_max=10**(numpy.log10(Dnu_star/beta0)/beta1)

	return nu_max

# The main function that generate a set of parameters used to generate Lorentzian profiles
# Assumptions: 
#     - The frequencies of l=0, l=2 and l=3 modes follow exactly the asymtptotic relation for the p mdoes
#	  - The frequencies of the l=1 modes follow exactly the asymptotitc relation for the mixed modes
#	  - Widths of l=0, 2, 3 modes are rescaled using the synthetic relation from Appourchaux+2014 applied to the solar profile
#	  - Widths of l=1 mixed modes are defined using the ksi function, scaled using l=0 modes. Thus l=0 modes fixes the upper
#		limit of the mode width
#	  - Heights of l=0 modes are rescaled using the measured heights of the solar modes
#	  - Heights of l=1 mixed modes are set equal to l=0 modes, modulo the bolometric visbility. 
#	    		WARNING: THIS IS AN INACCURATE ASSUMPTION ACCORDING TO MY PAPER. BUT NO FORMAL EXPRESSION 
#			                        IS KNOWN SO FAR SO THAT THIS IS THE BEST I CAN DO
#	  - Bolometric visibilities in height are assumed to be 1, 1.5, 0.5, 0.07 for l=0,1,2,3 respectively
#	  - Splitting is not implemented at all (zero-rotation)
# 	  Warning for widths and heights: if fmin/fmax is too large, you may have an error as an extrapolation 
#									  would be required, which I forbid by code. An fmin/fmax that englobes
#									  10 radial orders should not pose any issue.
# Input: 
#	Dnu_star: Large separation for the p modes
#   epsilon_star: phase offset term for the asymtptotic relation of the p modes
#   D0_star: D0 term, sensitive to core properties and in the asymtptoic relation for the p modes
#   DP1_star: The period spacing for l=1 g modes 
#	alpha_star: The phase offset term for the asymptotic relation of the g modes
#   q_star: Coupling coeficient between p and g modes
#	fmin: Minimum frequency for the modes that should be included in the calculations
#   fmax: Maximum frequency for the modes that should be included in the calculations
# Outputs:
#	nu_lx: Frequencies of the l=x modes. x is between 0 and 3
#	nu_p_l1: Base p modes frequencies used to build the frequencies for the l=1 mixed modes
#	nu_g_l1: Base p modes frequencies used to build the frequencies for the l=1 mixed modes
#	width_lx: Widths of the l=x modes. x is between 0 and 3
#   height_lx: Heights of the l=x modes. x is between 0 and 3 
def make_synthetic_asymptotic_star(Teff_star, numax_star, Dnu_star, epsilon_star, D0_star, DP1_star, alpha_star, q_star, fmin, fmax, Hmax_l0=1., Gamma_max_l0=1., rot_env_input=-1):

	# Fix the resolution to 4 years (converted into microHz)
	resol=1e6/(4*365.*86400.) 

	#Fix the bolometric visibilities
	Vl=[1, 1.5, 0.5, 0.07]

	# ----- l=0 modes -----
	# This section generate l=0 modes following the asymptotic relation of p modes, and make
	# rescaled width and height profiles for the star using the solar width and height profiles
	nmin=int(fmin/Dnu_star - epsilon_star)
	nmax=int(fmax/Dnu_star + epsilon_star)
	nu_l0=[]
	for en in range(nmin, nmax):
		tmp=solver_mm.asympt_nu_p(Dnu_star, en, epsilon_star, 0, D0_star)
		nu_l0.append(tmp)
	nu_l0=numpy.asarray(nu_l0)
	#width0, height0=width_height_sun_rescaled(nu_l0, Dnu_star, n_numax)
	width_l0, height_l0=width_height_MS_sun_rescaled(nu_l0, Dnu_star, numax_star)
	height_l0=height_l0*Hmax_l0

	width_l0=width_l0*Gamma_max_l0
	
	# Define an interpolation base to be used when generating maximum widths of l=1 mixed modes
	# and for heights and widths of l=2,3 modes	
	int_fct_h0 = interpolate.interp1d(nu_l0, height_l0)
	int_fct_w0 = interpolate.interp1d(nu_l0, width_l0)
	# ------------------------

	# ------- l=1 modes ------
	# Use the solver to get mixed modes
	el=1
	nu_m_l1, nu_p_l1, nu_g_l1=solver_mm.solve_mm_asymptotic(Dnu_star, epsilon_star, el, D0_star, DP1_star, alpha_star, q_star, fmin, fmax, resol, returns_pg_freqs=True)
	# Filter solutions that endup at frequencies higher/lower than the nu_l0 because we will need to extrapolate height/widths otherwise...
	posOK=numpy.where(numpy.logical_and(nu_m_l1 >= min(nu_l0), nu_m_l1 <= max(nu_l0)))
	nu_m_l1=nu_m_l1[posOK]

	# Generating widths profiles for l=1 modes using the ksi function
	Dnu_p=numpy.repeat(Dnu_star, len(nu_p_l1))
	DPl=numpy.repeat(DP1_star, len(nu_g_l1))

	ksi_pg=ksi_fct2(nu_m_l1, nu_p_l1, nu_g_l1, Dnu_p, DPl, q_star) # assunme Dnu_p, DPl and q constant
	h1_h0_ratio=h_l_rgb(ksi_pg) # WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019)
#	print('Len(h1_h0_ratio):',len(h1_h0_ratio))

	height_l1p=int_fct_h0(nu_m_l1)
	height_l1p=height_l1p*Vl[1] 
#	print('Len(Height_l1p):', len(height_l1p))
	height_l1=h1_h0_ratio * height_l1p 
	#width_l1=gamma_l_fct1(nu_m_l1, nu_p_l1, nu_g_l1, Dnu_p, DPl, q_star, nu_l0, width_l0, el, hl_h0_ratio=1.)
	width_l1=gamma_l_fct2(ksi_pg, nu_m_l1, nu_l0, width_l0, h1_h0_ratio, el)

	# Generating splittings with a two averaged rotation rates
	if rot_env_input <=0:
		# Determine the envelope rotation
		rot_env_input=rot_envelope()

	rot_env, rot_c, rot_env_true, rot_c_true=rot_2zones_v1(rot_env_input, numax_star, Teff_star, sigma_logg=0.1)
	a1_l1=dnu_rot_2zones(ksi_pg, rot_env, rot_c)

	#print(rot_env, rot_c, rot_env_true, rot_c_true)
	#print('a1_l1:')
	#print(a1_l1)
	#print('size(a1_l1)', len(a1_l1))
	#exit()
	# ------- l=2 modes -----
	nu_l2=[]
	for en in range(nmin, nmax):
		tmp=solver_mm.asympt_nu_p(Dnu_star, en, epsilon_star, 2, D0_star)
		nu_l2.append(tmp)
	nu_l2=numpy.asarray(nu_l2)

	# Filter solutions that endup at frequencies higher/lower than the nu_l0 because we will need to extrapolate height/widths otherwise...
	posOK=numpy.where(numpy.logical_and(nu_l2 >= min(nu_l0), nu_l2 <= max(nu_l0)))
	nu_l2=nu_l2[posOK]

	height_l2=int_fct_h0(nu_l2)
	height_l2=height_l2*Vl[2]
	width_l2=int_fct_w0(nu_l2)

	# Assume that the l=2 modes are only sensitive to the envelope rotation
	a1_l2=numpy.repeat(rot_env, len(nu_l2))

	# ------ l=3 modes ----
	nu_l3=[]
	for en in range(nmin, nmax):
		tmp=solver_mm.asympt_nu_p(Dnu_star, en, epsilon_star, 3, D0_star)
		nu_l3.append(tmp)
	nu_l3=numpy.asarray(nu_l3)

	posOK=numpy.where(numpy.logical_and(nu_l3 >= min(nu_l0), nu_l3 <= max(nu_l0)))
	nu_l3=nu_l3[posOK]

	height_l3=int_fct_h0(nu_l3)
	height_l3=height_l3*Vl[3]
	width_l3=int_fct_w0(nu_l3)

	# Assume that the l=2 modes are only sensitive to the envelope rotation
	a1_l3=numpy.repeat(rot_env, len(nu_l3))

	return nu_l0, nu_p_l1, nu_g_l1, nu_m_l1, nu_l2, nu_l3, width_l0, width_l1, width_l2, width_l3, height_l0, height_l1, height_l2, height_l3, a1_l1, a1_l2, a1_l3


# Read a basic configuration file that contains the setup to create the whole set
# of parameters for a synthetic spectra. Those parameters are returned in an output file.
def main_star_generator(config_file='star_params.global', output_file='star_params.modes', output_file_range='star_params.range'):

	comments=[]
	setup=[]
	f=open(config_file, 'r')
	for line in f:
		s=line.split()
		if s[0] != '#':
			setup.append(s)
		else:
			comments.append(line)
	f.close()

	setup_pmodes=numpy.asarray(setup[0], dtype='double')
	setup_gmodes=numpy.asarray(setup[1], dtype='double')
	coupling=numpy.double(setup[2][0])
	Ncoef=numpy.double(setup[2][1]) # Defines the number of modes on each side of numax
	Hmax_l0=numpy.double(setup[2][2]) # Maximum height for the l=0. Heights of other modes are set using visibilities (make_synthetic_asymptotic_star()).
	Gamma_max_l0=numpy.double(setup[2][3]) # Width at numax for l=0
	Dnu_star=numpy.double(setup_pmodes[0])
	numax_star=numax_from_stello2009(setup_pmodes[0])
	fmin=numax_star -Ncoef*Dnu_star
	fmax=numax_star +Ncoef*Dnu_star

	write_range_modes(numax_star, Dnu_star, Ncoef, output_file=output_file_range)
		
	nu_l0, nu_p_l1, nu_g_l1, nu_m_l1, nu_l2, nu_l3, width_l0, width_m_l1, width_l2, width_l3, height_l0, height_m_l1, height_l2, height_l3, a1_l1, a1_l2, a1_l3=make_synthetic_asymptotic_star(Teff_star, numax_star, setup_pmodes[0], setup_pmodes[1], setup_pmodes[2], setup_gmodes[0], setup_gmodes[1], coupling, fmin, fmax, Hmax_l0=Hmax_l0, Gamma_max_l0=Gamma_max_l0)
	fout=open(output_file, 'w')
	fout.write('# Configuration of mode parameters. This file was generated by bump_DP.py (external python program)')
	fout.write('# Input mode parameters. degree / freq / H / W / splitting a1 / eta / a3 /  b (mag) / alfa (mag) / asymetry / inclination')

	el=0
	for en in range(len(nu_l0)):
		nu=nu_l0[en]
		w=width_l0[en]
		h=height_l0[en]
		a1=0
		eta=0
		a3=0
		asym=0
		inc=0
		b=0
		alfa=0
		stri='{0}   {1:.5f}   {2:.5f}   {3:.5f}   {4:.2f}   {5:.10f}    {6:0.5f}    {7}    {8}   {9}    {10}'.format(el, nu, h, w, a1, eta, a3, b, alfa, asym, inc)
		fout.write(stri)
		fout.write('\n')
	el=1
	for en in range(len(nu_m_l1)):
		nu=nu_m_l1[en]
		w=width_m_l1[en]
		h=height_m_l1[en]
		a1=a1_l1[en]
		eta=0
		a3=0
		asym=0
		inc=0
		b=0
		alfa=0
		stri='{0}   {1:.5f}   {2:.5f}   {3:.5f}   {4:.2f}   {5:.10f}    {6:0.5f}    {7}    {8}   {9}    {10}'.format(el, nu, h, w, a1, eta, a3, b, alfa, asym, inc)
		fout.write(stri)
		fout.write('\n')
	el=2
	for en in range(len(nu_l2)):
		nu=nu_l2[en]
		w=width_l2[en]
		h=height_l2[en]
		a1=a1_l2[en]
		eta=0
		a3=0
		asym=0
		inc=0
		b=0
		alfa=0
		stri='{0}   {1:.5f}   {2:.5f}   {3:.5f}   {4:.2f}   {5:.10f}    {6:0.5f}    {7}    {8}   {9}    {10}'.format(el, nu, h, w, a1, eta, a3, b, alfa, asym, inc)
		fout.write(stri)
		fout.write('\n')
	el=3
	for en in range(len(nu_l3)):
		nu=nu_l3[en]
		w=width_l3[en]
		h=height_l3[en]
		a1=a1_l3[en]
		eta=0
		a3=0
		asym=0
		inc=0
		b=0
		alfa=0
		stri='{0}   {1:.5f}   {2:.5f}   {3:.5f}   {4:.2f}   {5:.10f}    {6:0.5f}    {7}    {8}   {9}    {10}'.format(el, nu, h, w, a1, eta, a3, b, alfa, asym, inc)
		fout.write(stri)
		fout.write('\n')
	fout.close()
	print('        ===> frequencies, heights, widths of l=0,1,2,3 p and mixed modes created successfully...')


def write_range_modes(numax, Dnu, Ncoef, output_file='star_params.range'):

	fmin=numax-Ncoef*Dnu
	fmax=numax+Ncoef*Dnu

	fout=open(output_file, 'w')
	fout.write('# min and max frequencies for the relevant modes. This file was generated by bump_DP.py (external python program)')
	fout.write('\n')
	stri='{0:.3f}   {1:.3f} '.format(fmin, fmax)
	fout.write(stri)
	fout.write('\n')
	fout.close()

# A function that allows you to test and visualise the results from make_synthetic_asymptotic_star()
def test_asymptotic_star(Dnu_star=20, DP1_star=80, q_star=0.15, Teff_star=5300.):	
	# Define global Pulsation parameters
	el=1.
	#Dnu_star=30. # RGB star
	#DP1_star=120. # RGB star
	#Dnu_star=55 #15 # microHz
	#DP1_star=350 #85 # microHz, typical for a RGB with Dnu_p=10 (see Fig. 1 of Mosser+2014, https://arxiv.org/pdf/1411.1082.pdf)
	#q_star=0.15 # Fix the coupling term

	# Parameters for p modes that follow exactly the asymptotic relation of p modes
	D0_star=Dnu_star/100. 
	epsilon_star=0.4

	# Parameters for g modes that follow exactly the asymptotic relation of g modes for a star with radiative core
	alpha_star=0.

	# Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	numax_star=numax_from_stello2009(Dnu_star)
	fmin=numax_star - 5*Dnu_star
	fmax=numax_star + 5*Dnu_star

	nu_l0, nu_p_l1, nu_g_l1, nu_m_l1, nu_l2, nu_l3, width_l0, width_m_l1, width_l2, width_l3, height_l0, height_l1, height_l2, height_l3, a1_l1, a1_l2, a1_l3=make_synthetic_asymptotic_star(Teff_star, numax_star, Dnu_star, epsilon_star, D0_star, DP1_star, alpha_star, q_star, fmin, fmax)

	plt.plot(nu_l0, width_l0, color='Black', linestyle='--')
	plt.plot(nu_l0, height_l0, color='Blue', linestyle='--')
	plt.plot(nu_m_l1, height_l1, color='Orange', linestyle='--')
	plt.plot(nu_l0, height_l0*1.5, color='Orange', linestyle='--')
	
	plt.plot(nu_m_l1, width_m_l1, color='Red', linestyle='--')
	plt.plot(nu_m_l1, width_m_l1, 'ro')
	plt.axvline(x=numax_star, color='red', linestyle='--')
	plt.axvline(x=numax_star - 5*Dnu_star, color='green', linestyle='--')
	plt.axvline(x=numax_star + 5*Dnu_star, color='green', linestyle='--')
	plt.show()

	# For testing ksi_fct2
	nu=numpy.linspace(fmin ,fmax, 10000)

	Dnu_p=numpy.repeat(Dnu_star, len(nu_p_l1))
	DPl=numpy.repeat(DP1_star, len(nu_g_l1))

	ksi=ksi_fct2(nu, nu_p_l1, nu_g_l1, Dnu_p, DPl, q_star)
	ksi_nu_m=ksi_fct2(nu_m_l1, nu_p_l1, nu_g_l1, Dnu_p, DPl, q_star, norm_method='slow')
	plt.plot(nu, ksi, color='Orange')
	plt.plot(nu_m_l1, ksi_nu_m, 'ro')
	for p in nu_p_l1:
		plt.axvline(x=p, color='blue', linestyle='--')
		print(p)
	plt.show()

	plt.plot(nu_m_l1, a1_l1, color='Red', linestyle='--')
	plt.plot(nu_m_l1, a1_l1, 'ro')
	#plt.plot(nu_m_l1, ksi_nu_m, 'ro')
	for p in nu_p_l1:
		plt.axvline(x=p, color='blue', linestyle='--')
		print(p)
	plt.show()


#main_star_generator(config_file='star_params.global', output_file='star_params.modes')


