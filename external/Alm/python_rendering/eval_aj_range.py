# eval_aj_range.py       : Contains a few functions that read an Alm python grid (can be created by make_Alm_grid.py) and use it give an idea of what should be the prior range for a given star as a function of some stellar parameters. The main function is 
import numpy as np
from activity import Qlm
from activity import eta0_fct
from acoefs import eval_acoefs
from termcolor import colored

def load_Alm_grids(dir_grid, gridfiles, lmax=2):
	Alm=[]
	if lmax>=1:
		Alm_grid_l1=np.load(dir_grid + gridfiles[0])
		Alm.append(Alm_grid_l1['Alm'])# Note: ftype = 'gauss' or 'gate' depends on the Alm grid file. No need to specify it
	if lmax >=2:
		Alm_grid_l2=np.load(dir_grid + gridfiles[1])
		Alm.append(Alm_grid_l2['Alm'])# Note: ftype = 'gauss' or 'gate' depends on the Alm grid file. No need to specify it
	if lmax >=3:
		Alm_grid_l3=np.load(dir_grid + gridfiles[2])
		Alm.append(Alm_grid_l3['Alm'])# Note: ftype = 'gauss' or 'gate' depends on the Alm grid file. No need to specify it	
	if lmax >=1:
		return Alm_grid_l1['theta'], Alm_grid_l1['delta'], Alm
	else:
		print('Warning: lmax<1. No grid can be retrieved')
		return [], [],[]

def numax_fct(Mass=1, Radius=1, Teff=5777.):
	numax_sun=3100.
	Teff_sun=5777.
	numax=numax_sun * Mass/(Radius**2 *np.sqrt(Teff/Teff_sun))
	return numax

def update_min(aj, aj_ref):
	if aj_ref > aj:
		return aj
	else:
		return aj_ref

def update_max(aj, aj_ref):
	if aj_ref < aj:
		return aj
	else:
		return aj_ref

#Assumptions: nu_max is derived from the requested Dnu_star parameter using the relation from Stello+2009. 
#	Dnu ~ 0.263*numax^0.77 (no uncertainty implemented here)
def numax_from_stello2009(Dnu_star):
	# Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	beta0=0.263; # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	beta1=0.77; # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	nu_max=np.power(10, np.log10(Dnu_star/beta0)/beta1)
	return nu_max

def Dnu_from_stello2009(numax):	
	# Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	beta0=0.263; # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	beta1=0.77; # according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	Dnu=beta0*np.power(numax, beta1)
	return Dnu

def eval_aj_mean_range(numax=None, Dnu=None, Mass=None, Radius=None, Teff=None, epsilon_nl=5*1e-4, a1=1000, dir_core='/Users/obenomar/tmp/test_a2AR/'):
	'''
	Evaluate the full range of possible of aj E [2,4] coefficients for a given star of Mass=Mass, Radius=Radius and Teff=Teff 
	and given an activity level epsilon_nl and a rotation a1
	If one wants a purely observational approach, it can provides directly numax (instead of M, R, Teff) and/or Dnu
	If only numax is provided, Dnu is determined using scaling relation with numax from Stello+2009. Same if only Dnu is provided
	This is usefull to set the range of simulations for aj and for setting priors on aj odds coefficient while fitting
	Note: This function **only** considers the mean value of aj over l and at numax
	'''
	R_sun=6.96342e5 #in km
	M_sun=1.98855e30 #in kg
	if numax != None and Dnu == None:
		print('Calculating Dnu(numax)...')
		Dnu=Dnu_from_stello2009(numax)
		print('Dnu = ', numax)
		eta0=eta0_fct(Dnu=Dnu)
	#
	if numax == None and Dnu != None:
		print('Calculating numax(Dnu)...')
		numax=numax_from_stello2009(Dnu)
		eta0=eta0_fct(Dnu=Dnu)
	#
	if numax != None and Dnu != None:
		print('Using numax and Dnu...')
		eta0=eta0_fct(Dnu=Dnu)
	#	
	if numax == None and Dnu == None and Mass != None and Radius != None and Teff !=None:
		print('Using M,R,Teff...')
		numax=numax_fct(Mass=Mass, Radius=Radius, Teff=Teff)
		Volume=np.pi*(Radius*R_sun)**3 * 4./3
		rho=Mass*M_sun/Volume * 1e-12 # Conversion into g/cm-3
		#print('rho =', rho)
		eta0=eta0_fct(rho=rho)
	if numax == None and Dnu == None and (Mass == None or Radius == None or Teff == None):
		print('Error in eval_aj_mean_range: You need to provide (M, R, Teff) or numax and/or Dnu')
		print('                             Please read the description of the function before using it')
		exit()
	
	print('numax =', numax)
	#print('eta0  =', eta0)
	lmax=3 # We do not consider a6... we evaluate the range for a2 and a4 only
	dir_grids=dir_core+"/grids/gate/0.25deg_resol/" # High precision grid
	gridfiles=['grid_Alm_1.npz', 'grid_Alm_2.npz', 'grid_Alm_3.npz'] # Fine grids
	#
	theta, delta, Alm= load_Alm_grids(dir_grids, gridfiles, lmax=lmax)
	#
	a2_min=99999
	a2_max=-99999
	a4_min=99999
	a4_max=-99999
	a6_min=99999
	a6_max=-999999
	for j in range(len(theta)):
		for k in range(len(delta)):
			aj=np.zeros(6)
			for l in range(1,lmax+1):
				nu_nlm=[]
				for m in range(-l, l+1):
					# perturvation from AR		
					dnu_AR=numax*epsilon_nl*Alm[l-1][m+l][j,k] # All Alm(theta, delta) are in [j,k]
					# perturbation from CF
					dnu_CF=eta0*numax * (a1*1e-9)**2 * Qlm(l,m)
					nu_nlm.append(numax + dnu_AR + dnu_CF)
					#print('(j,k) = (', j, ',', k,')  theta =', theta[j]*180./np.pi, "  delta =", delta[k]*180./np.pi, " dnu_AR =", dnu_AR*1000)
				a=eval_acoefs(l, nu_nlm)
				# Averaging over l
				for o in range(len(aj)):
					if o < 2:
						aj[o]=aj[o] + a[o]/lmax
					if o >=2 and o < 4:
						if l>=2:
							aj[o]=aj[o] + a[o]/(lmax-1)
					if o >=3 and o <6:
						if l>=3:
							aj[o]=aj[o] + a[o]/(lmax-2)
#					if j==10 and k==10: # For debug only
#						print('a[{}] = {}'.format(o, a[o]))
#			if j==10 and k==10 and l==lmax: # For debug only
#				print('aj = ', aj)
#				exit()
			a2_min=update_min(aj[1], a2_min)
			a2_max=update_max(aj[1], a2_max)
			a4_min=update_min(aj[3], a4_min)
			a4_max=update_max(aj[3], a4_max)
			a6_min=update_min(aj[5], a6_min)
			a6_max=update_max(aj[5], a6_max)
	return [a2_min*1e3, a2_max*1e3], [a4_min*1e3, a4_max*1e3], [a6_min*1e3, a6_max*1e3] # aj converted into nHz


def aj_mean_range_v1():
	'''
	This function is to explore specific range for some stars. It can for example be used to evaluate 
	The prior that has to be used before fitting a star for which M, R, Teff is known.
	'''
	R_sun=6.96342e5 #in km
	M_sun=1.98855e30 #in kg
	V=np.pi*R_sun**3 * 4./3
	rho_sun=M_sun/V * 1e-12 # Conversion into g/cm-3
	# Table of M, R derived from Bellinger+2019 (https://www.aanda.org/articles/aa/pdf/2019/02/aa34461-18.pdf)
	# Only Sun, Min and Max of his range are kept
	# Teff are comming from Pinsonault+2012 or Davies+2016 (or Campante+2015)
	
	KIC=['12069424']
	M=[1.056]
	R=[1.213]
	T=[5825]
	rho=[1.531]
	a1=[450]
	epsilon_nl=5e-4 # Solar Value
	print('epsilon_nl :', epsilon_nl)
	
	''''
	KIC=['12069449']
	M=[1.07]
	R=[1.127]
	T=[5750]
	#rho=[1.531]
	a1=[500]
	epsilon_nl=5e-4 # Solar Value
	'''
	for j in range(0,1):
		print('j=', j)
		print(colored('KIC =', 'red'), KIC[j])
		print('Mass: ', M[j],  '   Radius: ', R[j],   '    Teff: ', T[j])
		for dnu in a1:
			print(colored('a1         :'+str(dnu), 'red'))
			a2_range, a4_range, a6_range=eval_aj_mean_range(Mass=M[j], Radius=R[j], Teff=T[j], epsilon_nl=epsilon_nl, a1=dnu)
			print(colored('a2_range:', 'red'), a2_range, ' (nHz) | ', colored(100*np.abs(a2_range[0])/dnu, 'blue'),  ' , ', colored(np.abs(100*a2_range[1])/dnu, 'blue'), '  (%)')
			print(colored('a4_range:', 'red'), a4_range, ' (nHz) | ', colored(100*np.abs(a4_range[0])/dnu, 'blue'),  ' , ', colored(np.abs(100*a4_range[1])/dnu, 'blue'), '  (%)')
			print(colored('a6_range:', 'red'), a6_range, ' (nHz) | ', colored(100*np.abs(a6_range[0])/dnu, 'blue'),  ' , ', colored(np.abs(100*a6_range[1])/dnu, 'blue'), '  (%)')
			print(colored('----', 'red'))

#aj_mean_range_v1()
