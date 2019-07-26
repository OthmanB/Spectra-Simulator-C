# This contains all the functions for solving the asymptotic relation
# of the mixed modes, as they have been tested during their development
# All this arise from reading few papers from Benoit Mosser and 
# The PhD thesis from Charlotte Gehand:
# https://arxiv.org/pdf/1203.0689.pdf
# https://arxiv.org/pdf/1411.1082.pdf
# https://tel.archives-ouvertes.fr/tel-02128409/document

# Examples and tests function have been built using asymptotic relations.
# But note that they should be applicable to ANY set of value of:
#	 nu_p(nu), nu_g(nu), Dnu_p(nu) and DPl(nu) (meaning handling glitches)
import numpy
#import matplotlib.pyplot as plt
from scipy import interpolate

def asympt_nu_p(Dnu_p, np, epsilon, l, D0):
	return (np + epsilon + l/2.)*Dnu_p - l*(l+1)*D0

def asympt_nu_g(DPl, ng, alpha):
	Pl=(ng + alpha)*DPl
	return 1e6/Pl

def pnu_fct(nu, nu_p):
	try:
		l=len(nu)
		pnu=nu - numpy.repeat(nu_p, len(nu))
	except:
		l=1
		pnu=nu - nu_p

	return pnu

def gnu_fct(nu, nu_g, Dnu_p, DPl, q):
	X=numpy.pi * (1./nu - 1./nu_g)*1e6 / DPl
	gnu=Dnu_p*numpy.arctan(q * numpy.tan(X)) / numpy.pi
	return gnu


def rem_doubles(arr=[1.,1.01,2.,2.2,2.4,2.41,2.43,2.45,3., 4., 5., 5.5, 5.51, 5.52, 5.53, 5.5301, 6.], tol=0.1):
	L=len(arr)
	sol=[]
	i=0
	endofarray=False
	while i < L-1:
		j0=i
		j=j0
		delta=0
		while delta <= tol and endofarray==False:
			if j == L-1:	
				endofarray=True
			else:
				delta=numpy.abs(arr[j] - arr[j+1])
			j=j+1
		if endofarray == False:
			if j0 != j-1:
				sol.append(numpy.mean(arr[j0:j-1]))
			else:
				sol.append(numpy.mean(arr[j0]))
		else:
			sol.append(numpy.mean(arr[j0:]))
		i=j
	if arr[L-1] -arr[L-2] >= tol and endofarray == False:
		sol.append(arr[L-1])
	return sol


'''
This the main function that solves the mixed mode asymptotic relation
which is of the type p(nu) = g(nu)
This solver specifically solve the case:
      nu - nu_p = Dnu*arctan(q tan(1/(nu DPl) - 1/(nu_g*DPl)))
      It tries to find the intersect between p(nu) and g(nu)
      using an interpolation and over a range of frequency such that nu is in [numin, numax]
Parameters:
	- Mandatory: 
	     nu_p (double) : frequency of a given p-mode (in microHz)
	     nu_g (double): frequency of a given g-mode (in microHz)
	     Dnu_p (double): Large separation for p modes (in microHz)
	     DP1 (double): Period spacing for g modes (in seconds)
	     q (double): Coupling term (no unit, should be between 0 and 1)
	- Optional:
		numin (double): Minimum frequency considered for the solution (in microHz)
		numax (double): Maximum frequency considered for the solution (in microHz)
		resol (double): Base resolution for the interpolated base function. The interpolation may miss solutions 
		       if this is set too low. Typically, the resolution parameter should be higher than the
		       spectral resolution of your spectrum so that all resolved modes should be found.
		       This is also used for creating the nu axis for visualisation (in microHz).
		factor (double): Define how fine will be the new tiny grid used for performing the interpolation. This is important
				to avoid extrapolation (which is forbiden and will result in crash of the code). Typically, the default
				value factor=0.05 can compute mixed modes for frequency down to 80microHz. Going below requires a smaller factor

		returns_axis: If True, returns nu, pnu and gnu (see optional reutrns below). Mainly for debug
Returns:
	nu_m: An array with all solutions that match p(nu) = g(nu)
	nu (optional): The frequency axis used as reference for finding the intersection
	pnu (optional): The curve for p(nu)
	gnu (optional): The curve g(nu)
'''
def solver_mm(nu_p, nu_g, Dnu_p, DPl, q, numin=1, numax=1500., resol=1, returns_axis=False, verbose=False, factor=0.05):

	# Generate a frequency axis that has a fixed resolution and that span from numin to numax
	nu=numpy.linspace(numin, numax, num=(numax-numin)/resol)

	# Function p(nu) describing the p modes
	pnu=pnu_fct(nu, nu_p)

	# Function g(nu) describing the g modes 
	gnu=gnu_fct(nu, nu_g, Dnu_p, DPl, q)

	#plt.plot(nu, pnu, color='black')
	#plt.plot(nu, gnu, 'r')
	
	# Find when p(nu) = g(nu) by looking for solution of p(nu) - g(nu) = 0
	#     Method 1: Direct Interpolation... Works only for single solutions ==> Not used here
	#int_fct = interpolate.interp1d(pnu - gnu, nu)
	#nu_m=int_fct(0)
	#     Method 2: (a) Find indices close to sign changes for p(nu) - g(nu)
	#               (b) Then perform an iterative interpolation in narrow ranges
	#                   near the approximate solutions. How narrow is the range is defined
	#					by the resolution parameter resol, which in this case can be view
	#					as the minimum precision.
	nu_m=[]
	idx = numpy.argwhere(numpy.diff(numpy.sign(pnu - gnu))).flatten()
	#print(numpy.diff(numpy.sign(pnu - gnu)))
	#print('Indexes:', idx)
	for ind in idx:
		# Define a small local range around each of the best solutions
		range_min=nu[ind] - 2*resol
		range_max=nu[ind] + 2*resol
		# Redefine nu, pnu and gnu for that local range
		nu_local=numpy.linspace(range_min, range_max, num=(range_max-range_min)/(resol*factor))
		pnu_local=pnu_fct(nu_local, nu_p)
		gnu_local=gnu_fct(nu_local, nu_g, Dnu_p, DPl, q)	

		#plt.plot(nu_local, pnu_local, color='blue')
		#plt.plot(nu_local, gnu_local, color='purple')
		#plt.show()
		#exit()

		# Perform the interpolation on the local range and append the solution to the nu_m list
		int_fct = interpolate.interp1d(pnu_local - gnu_local, nu_local)
		nu_m_proposed=int_fct(0)
		ysol_gnu=gnu_fct(nu_m_proposed, nu_g, Dnu_p, DPl, q)
		ysol_pnu=pnu_fct(nu_m_proposed, nu_p)
		ratio=ysol_gnu/ysol_pnu
		if verbose == True:
			print('-------')
			print('nu_m:', nu_m_proposed)
			print('Ratio:', ratio)

		# Sometimes, the interpolator mess up due to the limits of validity for the atan function
		# The way to keep real intersection is to verify after interpolation that we really
		# have p(nu_m_proposed) = g(nu_m_proposed). We then only keeps solutions that satisfy
		# a precision criteria of 1%.
		if (ratio >= 0.999) and (ratio <= 1.001):
			nu_m.append(nu_m_proposed)

		#plt.plot(nu_local, pnu_local, color='b')
		#plt.plot(nu_local,gnu_local, color='r')
		#plt.plot(int_fct(0), ysol_gnu, 'ro')
		#plt.plot(int_fct(0), ysol_pnu, 'bx')
		#plt.show()
	
	nu_m=numpy.asarray(nu_m)
	
	ysol_gnu=gnu_fct(nu_m, nu_g, Dnu_p, DPl, q)
	#ysol_pnu=pnu_fct(nu_m, nu_p)
	ysol=ysol_gnu

	#print('nu_m:' , nu_m)

	if returns_axis == True:
		return nu_m, ysol, nu,pnu, gnu
	else:
		return nu_m, ysol

# This function uses solver_mm to find solutions from a spectrum
# of pure p modes and pure g modes following the asymptotic relations (no glitch, no second order)
def solve_mm_asymptotic(Dnu_p, epsilon, el, D0, DPl, alpha, q, fmin, fmax, resol, returns_pg_freqs=True, verbose=False):

	returns_axis=True

	# Use fmin and fmax to define the number of pure p modes and pure g modes to considered
	np_min=int(numpy.floor((fmin + D0)/Dnu_p - epsilon - el/2))
	np_max=int(numpy.ceil((fmax + D0)/Dnu_p - epsilon - el/2))

	ng_min=int(numpy.floor(1e6/(fmax*DPl) - alpha))
	ng_max=int(numpy.ceil(1e6/(fmin*DPl) - alpha))

	if fmin <= 100:
		fact=0.01
	else:
		fact=0.04
	nu_p_all=[]
	nu_g_all=[]
	nu_m_all=[]
	for np in range(np_min, np_max):
		for ng in range(ng_min, ng_max):
			nu_p=asympt_nu_p(Dnu_p, np, epsilon, el, D0)
			nu_g=asympt_nu_g(DPl, ng, alpha)
			nu_p_all.append(nu_p)
			nu_g_all.append(nu_g)
			nu_m, ysol, nu,pnu, gnu=solver_mm(nu_p, nu_g, Dnu_p, DPl,  q, numin=nu_p - Dnu_p, numax=nu_p + Dnu_p, resol=resol, returns_axis=returns_axis, factor=fact)
			if verbose == True:
				print('==========================================')
				print('nu_p: ', nu_p)
				print('nu_g: ', nu_g)
				print('solutions nu_m: ', nu_m)
			for sol in nu_m:
				nu_m_all.append(sol)

	nu_p_all=numpy.asarray(nu_p_all)
	nu_g_all=numpy.asarray(nu_g_all)

	#print('Before removing duplicates: ', nu_m_all)

	# Cleaning doubles: Solution 1 assuming exact matches (Fist Filter)
	result=[]
	for sol in nu_m_all:
		if sol not in result:
			result.append(sol)
	result=numpy.asarray(result)
	#print('After removing duplicates: ', result)

	# Cleaning doubles: Solution 2 removing using a tolerance range (Second Filter)
	result=numpy.sort(result)
	nu_m_final=rem_doubles(arr=result, tol=0.01)
	nu_m_final=numpy.asarray(nu_m_final)

	rp=[]
	rg=[]
	if returns_pg_freqs == True:
		for sol in nu_p_all:
			if sol not in rp:
				rp.append(sol)
		for sol in nu_g_all:
			if sol not in rg:
				rg.append(sol)

		rp=numpy.asarray(rp)
		rg=numpy.asarray(rg)

		return nu_m_final, rp, rg
	else:
		return nu_m_final

# Function to test solver_mm()
# This is a typical RGB case, with  density of g modes >> density of p modes
def test_rgb_solver_mm():

	Dnu_p=15 # microHz
	DP1= 80 # microHz, typical for a RGB with Dnu_p=10 (see Fig. 1 of Mosser+2014, https://arxiv.org/pdf/1411.1082.pdf)

	# Generate a p-mode that follow exactly the asymptotic relation of p modes
	D0=Dnu_p/100. 
	epsilon=0.4
	np=10.
	nu_p=(np + epsilon + 1./2.)*Dnu_p - 2*D0
	# Generate a g-mode that follow exactly the asymptotic relation of g modes for a star with radiative core
	ng=50
	alpha=0.01
	nu_g=1e6/(ng*DP1)

	# Use the solver
	q=0.1 # Fix the coupling term
	nu_m, ysol, nu,pnu, gnu=solver_mm(nu_p, nu_g, Dnu_p, DP1, q, numin=nu_p-Dnu_p/2, numax=nu_p + Dnu_p/2, resol=0.01, returns_axis=True)
	#nu_m, ysol, nu,pnu, gnu=solver_mm(nu_p, nu_g, Dnu_p, DP1, q, numin=1, numax=1000, resol=0.01, returns_pg_freqs=True)

	# Plot the outputs to check intersection are properly found
	plt.plot(nu, pnu, color='b')
	plt.plot(nu,gnu, color='r')
	plt.plot(nu_m, ysol, 'ro')
	plt.show()

	return nu_m, ysol, nu,pnu, gnu
	

# Function to test solver_mm()
# This is a typical SG case, with  density of g modes << density of p modes
def test_sg_solver_mm():

	Dnu_p=60 # microHz
	DP1= 400 # microHz, typical for a RGB with Dnu_p=10 (see Fig. 1 of Mosser+2014, https://arxiv.org/pdf/1411.1082.pdf)

	# Generate a p-mode that follow exactly the asymptotic relation of p modes
	D0=Dnu_p/100. 
	epsilon=0.4
	np=10.
	nu_p=(np + epsilon + 1./2.)*Dnu_p - 2*D0
	# Generate a g-mode that follow exactly the asymptotic relation of g modes for a star with radiative core
	ng=10
	alpha=0
	nu_g=1e6/(ng*DP1)

	# Use the solver
	q=0.2 # Fix the coupling term
	nu_m, ysol, nu,pnu, gnu=solver_mm(nu_p, nu_g, Dnu_p, DP1, q, numin=nu_p - Dnu_p/2, numax=nu_p + Dnu_p/2, resol=0.01, returns_axis=True)

	# Plot the outputs to check that intersection are properly found
	plt.plot(nu, pnu, color='b')
	plt.plot(nu,gnu, color='r')
	plt.plot(nu_m, ysol, 'ro')
	plt.show()

	return nu_m, ysol, nu, pnu, gnu

# Function to test solve_mm_asymptotic
# The parameters are typical for a RGB in the g mode asymptotic regime
def test_rgb_asymptotic():

	# Define global Pulsation parameters
	el=1.
	Dnu_p=15.
	DPl=80.
	q=0.1 # Fix the coupling term

	# Parameters for p modes that follow exactly the asymptotic relation of p modes
	D0=Dnu_p/100. 
	epsilon=0.4

	# Parameters for g modes that follow exactly the asymptotic relation of g modes for a star with radiative core
	#ng=10
	alpha=0.01

	# Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	beta0=0.263 # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	beta1=0.77 # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	nu_max=10**(numpy.log10(Dnu_p/beta0)/beta1)

	fmin=nu_max - 4*Dnu_p
	fmax=nu_max + 4*Dnu_p

	# Fix the resolution to 4 years (converted into microHz)
	resol=1e6/(4*365.*86400.) 

	# Use the solver
	freqs_mixed=solve_mm_asymptotic(Dnu_p, epsilon, el, D0, DPl, alpha, q, fmin, fmax, resol, returns_pg_freqs=False, verbose=False)

	# ---- Making an echelle diagram ------
	np_min=int(numpy.floor((fmin + D0)/Dnu_p - epsilon - el/2))
	np_max=int(numpy.ceil((fmax + D0)/Dnu_p - epsilon - el/2))

	ng_min=int(numpy.floor(1e6/(fmax*DPl) - alpha))
	ng_max=int(numpy.ceil(1e6/(fmin*DPl) - alpha))

	np=numpy.linspace(np_min, np_max, np_max - np_min + 1)
	freqs_l0=asympt_nu_p(Dnu_p, np, epsilon, 0, D0)
	freqs_l1_p=asympt_nu_p(Dnu_p, np, epsilon, 1, D0)

	# ---- This is if you want an echelle diagram in Frequency ----
	xaxis=numpy.linspace(0, DPl, 5)
	ng=numpy.linspace(ng_min, ng_max, ng_max - ng_min + 1)
	freqs_l1_g=asympt_nu_g(DPl, ng, alpha)

	plt.plot(1e6/freqs_l0 % DPl, 1e6/freqs_l0, 'bx')
	#plt.plot(1e6/freqs_l1_p % DPl, 1e6/freqs_l1_p, color='blue')
	plt.plot(1e6/freqs_mixed % DPl, 1e6/freqs_mixed, 'ro')
	plt.plot(1e6/freqs_l1_g % DPl, 1e6/freqs_l1_g, 'bo')
	for fp in freqs_l1_p:
		plt.plot(xaxis, 1e6/fp - xaxis, 'b')

	plt.xlim((-1,DPl*1.06))
	# --------------------------------------------------------------

	print(' --- Lenghts ----')
	print('L(nu_g): ', len(freqs_l1_g))
	print('L(nu_p): ', len(freqs_l1_p))
	print('L(nu_m): ', len(freqs_mixed))

	plt.show()

	return freqs_mixed

# Function to test solve_mm_asymptotic
# The parameters are typical for a RGB in the g mode asymptotic regime
def test_sg_asymptotic():

	# Define global Pulsation parameters
	el=1.
	Dnu_p=60 # microHz
	DPl= 400 # microHz, typical for a RGB with Dnu_p=10 (see Fig. 1 of Mosser+2014, https://arxiv.org/pdf/1411.1082.pdf)
	q=0.1 # Fix the coupling term

	# Parameters for p modes that follow exactly the asymptotic relation of p modes
	D0=Dnu_p/100. 
	epsilon=0.4

	# Parameters for g modes that follow exactly the asymptotic relation of g modes for a star with radiative core
	#ng=10
	alpha=0.

	# Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	beta0=0.263 # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	beta1=0.77 # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	nu_max=10**(numpy.log10(Dnu_p/beta0)/beta1)

	fmin=nu_max - 6*Dnu_p
	fmax=nu_max + 4*Dnu_p

	# Fix the resolution to 4 years (converted into microHz)
	resol=1e6/(4*365.*86400.) 

	# Use the solver
	freqs_mixed=solve_mm_asymptotic(Dnu_p, epsilon, el, D0, DPl, alpha, q, fmin, fmax, resol, returns_pg_freqs=False)


	# ---- Making an echelle diagram ------
	np_min=int(numpy.floor((fmin + D0)/Dnu_p - epsilon - el/2))
	np_max=int(numpy.ceil((fmax + D0)/Dnu_p - epsilon - el/2))

	ng_min=int(numpy.floor(1e6/(fmax*DPl) - alpha))
	ng_max=int(numpy.ceil(1e6/(fmin*DPl) - alpha))

	np=numpy.linspace(np_min, np_max, np_max - np_min + 1)
	freqs_l0=asympt_nu_p(Dnu_p, np, epsilon, 0, D0)
	freqs_l1_p=asympt_nu_p(Dnu_p, np, epsilon, 1, D0)


	# ---- This is if you want an echelle diagram in Frequency ----
	xaxis=numpy.linspace(0, Dnu_p, 5)
	ng=numpy.linspace(ng_min, ng_max, ng_max - ng_min + 1)
	freqs_l1_g=asympt_nu_g(DPl, ng, alpha)

	plt.plot(freqs_l0 % Dnu_p, freqs_l0, marker='^', color='black')
	plt.plot(freqs_l1_p % Dnu_p, freqs_l1_p, marker='s', color='blue')
	plt.plot(freqs_mixed % Dnu_p, freqs_mixed, 'ro')
	plt.plot(freqs_l1_g % Dnu_p, freqs_l1_g, 'bx')
	for fg in freqs_l1_g:
		plt.plot(xaxis, fg - xaxis, 'r')
	# --------------------------------------------------------------

	print(' --- Lenghts ----')
	print('L(nu_g): ', len(freqs_l1_g))
	print('L(nu_p): ', len(freqs_l1_p))
	print('L(nu_m): ', len(freqs_mixed))

	plt.show()

	return freqs_mixed
