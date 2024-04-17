import numpy as np
import matplotlib.pyplot as plt
from scipy import special, integrate
from subprocess import Popen, PIPE
from acoefs import eval_acoefs

def Qlm(l,m):
	Dnl=2./3
	Qlm=(l*(l+1) - 3*m**2)/((2*l - 1)*(2*l + 3))
	return Qlm*Dnl

def eta0_fct(Dnu=None, rho=None, verbose=False):
	if rho == None:
		if Dnu !=None:
			Dnu_sun=135.1
			#numax_sun=3150.
			R_sun=6.96342e5 #in km
			M_sun=1.98855e30 #in kg
			rho_sun=M_sun*1e3/(4*np.pi*(R_sun*1e5)**3/3) #in g.cm-3
			rho=(Dnu/Dnu_sun)**2 * rho_sun
			if verbose == True:
				print('Computing rho using Dnu...')
				print(' rho =', rho)
		else:
			print('Error in eta0(): You need to provide at least one of the following argument:')
			print('                 - Dnu')
			print('                 - rho')
			exit()
	G=6.667e-8 # cm3.g-1.s-2
	#eta0=3./(4.*np.pi*rho*G) # second^2 # WRONG BECAUSE WE ASSUMED OMEGA ~ a1. IT SHOULD BE OMEGA ~ 2.pi.a1
	eta0=3.*np.pi/(rho*G)
	return eta0

def triangle_filter(theta, theta0, delta):
	'''
	A triangular filter which peaks at theta = theta0 and F=1
	Defined between theta = [0, pi/2]
	But I allow theta = [pi/2, pi] as well, TO BE USED FOR VISUALISATION ONLY
	For this model to make sense, note that theta0 must be in [0,pi/2]
	'''
	if np.min(theta) < 0:
		print('Error: Some values of theta are below 0. I do not allow this here for safety reason')
		exit()
	#
	try:
		l=len(theta)
	except:
		l=1
		theta=np.zeros(1, dtype=float) + theta
	F=np.zeros(l, dtype=float)
	# Lower part: Between 0<theta<theta0
	g1=np.where(np.bitwise_and(theta <= theta0, theta >=0))[0]
	a=2/delta # The slope
	b=1- a*theta0
	for g in g1:
		if (theta[g] - (theta0-delta/2)) > 0:
			F[g]=a*theta[g] + b
	# Higher part: Between theta0<theta<np.pi/2
	g2=np.where(np.bitwise_and(theta > theta0, theta <=np.pi/2))[0]
	a=-2/delta
	b=1-a*theta0
	for g in g2:
		if (theta[g]- (theta0 + delta/2)) < 0:
			F[g]=a*theta[g] + b
	#
	if np.max(theta) > np.pi/2:
		print('Warning: You use the function in the range above pi/2')
		print('         Suitable only for vidsualisation or debug purpose')
		print(theta)
		print('theta0 =', theta0)
		print('delta = ', delta)
	# Lower part: Between 0<theta<theta0
		g1bis=np.where(np.bitwise_and(theta <= np.pi - theta0, theta >=np.pi/2))[0]
		a=2/delta # The slope
		b=1- a*(np.pi - theta0)
		for g in g1bis:
			if (theta[g] - (np.pi - theta0-delta/2)) > 0:
				F[g]=a*theta[g] + b
	# Higher part: Between theta0<theta<np.pi/2
		g2=np.where(np.bitwise_and(theta > np.pi - theta0, theta <=np.pi))[0]
		a=-2/delta
		b=1-a*(np.pi - theta0 )
		for g in g2:
			if (theta[g]- (np.pi - theta0 + delta/2)) < 0:
				F[g]=a*theta[g] + b
	return F


def gauss_filter(theta, theta0, delta):
	F=np.exp(-(theta - theta0)**2/(2*delta**2)) + np.exp(-(theta - np.pi + theta0)**2/(2*delta**2)) 
	return F


def gate_filter(theta, theta0, delta):
	try:
		l=len(theta)
		F=np.zeros(l, dtype=float)
		g1=np.where(np.bitwise_and(theta >= (theta0 - delta/2), theta <= (theta0 + delta/2)))
		g2=np.where(np.bitwise_and(theta >= (np.pi - theta0 - delta/2), theta <= (np.pi - theta0 + delta/2)))
		for g in g1:
			F[g]=1
		for g in g2:
			F[g]=1
	except:
		if (theta >= (theta0 - delta/2)) and (theta <= (theta0 + delta/2)) or (theta >= (np.pi - theta0 - delta/2) and  theta <= (np.pi - theta0 + delta/2)):
			F=1
		else:
			F=0
	return F


def gauss_filter_cte(theta0, delta):
	'''
		A function that calculate what is the max value of a double-gaussian filter
		that is symetrical towards pi/2 and bounded between [0, pi] and of width delta
                     -                      -
		          -    -                  -   -
		         -       -              -      -
		        -         -            -        -
		       -            -         -          - 
		      -              -       -            -
		    -                  -   -                -
		  -                      -                    -
		--+----------+-----------+---------------------+----> theta
		  0         theta0      pi/2    pi-theta0      pi
	'''
	theta= np.linspace(0, np.pi, 100)
	F=gauss_filter(theta, theta0, delta)
	#F=F/max(F)
	return max(F)

def show_filters(theta0=np.pi/6, delta=np.pi/10):
	theta=np.linspace(0, np.pi, 1000)
	F1=gate_filter(theta, theta0, delta)
	F2=gauss_filter(theta, theta0, delta)
	F2=F2/gauss_filter_cte(theta0, delta)
	F3=triangle_filter(theta, theta0, delta)
	plt.plot(theta/np.pi, F3)
	plt.plot(theta/np.pi, F1)
	plt.plot(theta/np.pi, F2)
	plt.axvline(x=0.5, linestyle='--', color='red')
	plt.axhline(y=1, linestyle='--', color='red')
	plt.axhline(y=0, linestyle='--', color='red')
	plt.axvline(x=theta0/np.pi, color='green', linestyle='-.')
	plt.axvline(x=1 - theta0/np.pi, color='green', linestyle='-.')
	plt.xlabel('co-latitude (in pi units)')
	plt.show()

# -----

def tests_gate(l=1,m=1):
	theta_min=np.asarray([0, 30, 60.])*np.pi/180.
	theta_max=np.asarray([30, 60, 90])*np.pi/180.
	
	I_0_pi=0 # The full integral calculated in segments as per defined in theta_min / theta_max
	I_0_pi_filt=0 # Same as I_0_pi but when using the filtering function
	I_0_pi_cpp=0 # For the result of the CPP calculation
	for i in range(len(theta_min)):
		theta0=(theta_min[i] + theta_max[i])/2
		delta=theta_max[i] - theta_min[i]
		I0=integrate_Alm_gate(l, m, theta0, delta, phi_range=[0, 2*np.pi])
		I1=integrate_Alm_gate(l, m, np.pi-theta0, delta, phi_range=[0, 2*np.pi])
		Ifilt=integrate_Alm(l, m, [0, 2*np.pi], [0, np.pi], theta0, delta, ftype='gate')
		el, em, Ifilt_cpp=Alm_cpp(l, theta0, delta, "gate", raw=False, use2pi=False)
		I_0_pi=I_0_pi + I0[0]
		I_0_pi_filt=I_0_pi_filt + Ifilt[0]		
		I_0_pi_cpp=I_0_pi_cpp+Ifilt_cpp[m-l]
		print(" (theta_min, theta_max) = ({} , {})".format(theta_min[i], theta_max[i]))
		print("                 theta0 = ", theta0 )
		print("                  delta = ", delta )		
		print("   1. Calculation with direct bundaries: ")
		print("      North hemishpere:", I0[0])
		print("      South hemishpere:", I1[0])
		print("                 Total:", I0[0] + I1[0])
		print("   2. Calculation with Filter function (Python code): ")
		print("      North hemishpere: Not available")
		print("      South hemishpere: Not available")
		print("                 Total:", Ifilt[0])
		print("   3. Calculation with Filter function (CPP code): ")
		print("      North hemishpere: Not available")
		print("      South hemishpere: Not available")
		print("                 Total:", Ifilt_cpp[m-l]) # We pick only the compared m
		print("------")
		#plt.plot(theta_axis, Y)
	print(" Integral over the full range 0 - Pi :")
	print("    1. With direct bundaries (sum of all previous totals) :", 2*I_0_pi)
	print("    2. With the filter (sum of all previous totals) :", I_0_pi_filt)
	print("    3. With the CPP code (sum of all previous totals) :", I_0_pi_cpp)
	print("    4. Direct bundaries: ", integrate_Alm_gate(l, m, np.pi/2, np.pi, phi_range=[0, 2*np.pi]) ) # This is a full window fit [0, Pi]
	print("    5. Filter function : ", integrate_Alm(l, m, [0, 2*np.pi],[0, np.pi], np.pi/4, np.pi/2, ftype='gate' ) ) # Note that the difference in range with direct bundaries is NORMAL
	el, em, Int=Alm_cpp(l, np.pi/4, np.pi/2, "gate", raw=False, use2pi=False)
	print("    6. With the CPP code : ", Int[m-l] ) 

def test_integrate_ylm2(l):
    phi_range = [0, 2.*np.pi]
    theta_range = [0, np.pi/4.]
    print("Ylm2(", l, "):")
    for m in range(-l, l+1):
        integral=integrate_ylm2(l, m, phi_range, theta_range)
        print("(l=" , l, ", m=", m, ") =" , "   ", integral[0])


def Alm_triangle(_theta, _phi, _l, _m, _theta0, _delta):
	Y=Ylm2(_theta, _phi, _l, _m)
	F=triangle_filter(_theta, _theta0, _delta)
	return Y*F


def Alm_gate(_theta, _phi, _l, _m, _theta0, _delta):
	Y=Ylm2(_theta, _phi, _l, _m)
	F=gate_filter(_theta, _theta0, _delta)
	return Y*F


def Alm_gauss(_theta, _phi, _l, _m, _theta0, _delta):
	Y=Ylm2(_theta, _phi, _l, _m)
	F=gauss_filter(_theta, _theta0, _delta)
	Fmax=gauss_filter_cte(_theta0, _delta)
	return Y*F/Fmax

def Ylm2(_theta, _phi, _l, _m):
    _s = special.sph_harm(_m, _l, _phi, _theta)
    return (_s.real**2 + _s.imag**2) * np.sin(_theta)

def integrate_ylm2(l, m, phi_range, theta_range):
    result = integrate.dblquad(Ylm2,
         phi_range[0], phi_range[1], theta_range[0], theta_range[1], args=(l, m,))
    return result

def integrate_Alm(l, m, phi_range, theta_range, theta0, delta, ftype='gate'):
	if delta != 0:
		if ftype == 'triangle':
			result = integrate.dblquad(Alm_triangle,
		     phi_range[0], phi_range[1], theta_range[0], theta_range[1], args=(l, m, theta0, delta,))
		if ftype == 'gate':
			result = integrate.dblquad(Alm_gate,
		     phi_range[0], phi_range[1], theta_range[0], theta_range[1], args=(l, m, theta0, delta,))
		if ftype == 'gauss':
			result = integrate.dblquad(Alm_gauss,
		     phi_range[0], phi_range[1], theta_range[0], theta_range[1], args=(l, m, theta0, delta,))
		if ftype != 'gate' and ftype != 'gauss' and ftype != 'triangle':
			print("Wrong filter type argument: ")
			print("Use only:")
			print("    ftype='gate' or  ftype='gauss' or ftype='triangle'")
			print("The program will exit now ")
			exit()
	else:
		result=[0,0] # When delta is 0, obviously the result is 0
	return result

def integrate_Alm_gate(l, m, theta0, delta, phi_range=[0, 2*np.pi]):
	if delta != 0:
			result = integrate.dblquad(Ylm2,
		     phi_range[0], phi_range[1], theta0-delta/2, theta0+delta/2, args=(l, m,))
	else:
		result=[0,0] # When delta is 0, obviously the result is 0
	return result

def Alm_cpp(l, theta0, delta, ftype, raw=False, Alm_path='cpp_prg/'):
	# Function updated on 15/04/2023: To accomodate the new option format 
	# of the Alm program
	try:
		print([Alm_path + "./Alm", "-d ", l,
		   	"--theta0", theta0, "--delta", delta, 
			"-f", ftype])
		process = Popen([Alm_path + "./Alm", "-d", str(l),
		   	"--theta0", str(theta0), "--delta", str(delta), 
			"-f", ftype], stdout=PIPE, stderr=PIPE)
		(output, err) = process.communicate()
		exit_code = process.wait()
		#print(output)
		output=output.decode("utf-8") 
		if raw == False:
			r=output.split('\n')
			#config=[]
			l=[]
			m=[]
			Alm=[]
			for line in r:
				line=line.strip()
				#print("line =", line)
				try:
					if line[0] != "#":
						s=line.split()
						l.append(int(s[0]))
						m.append(float(s[1]))
						Alm.append(float(s[2]))
				except:
					err=True
			return np.asarray(l),np.asarray(m),np.asarray(Alm)
		else:
			return output, err
	except: # Handling processes that does not exist, aka, when  the Alm file is not available
		error=True
		print("Error: Could not execute the Alm C++ program. The most likely explanation is that it is not in the current directory")
		return -1, error

def Alm(l,m, theta0=np.pi/2, delta=2*8.4*np.pi/180, ftype='gate'):
    phi_range = [0, 2.*np.pi]
    theta_range = [0, np.pi/2]
    integral=integrate_Alm(l, m, phi_range, theta_range, theta0, delta, ftype=ftype)
    return integral[0]*2 # The final result is given over [0, pi], which is twice [0, np.pi/2] by definition (Alm axi-symetrical toward the Equator)


# Compute the a-coefficients for the theoretical model and provided key parameters of that model
# Use Alm_cpp instead of Alm in python... much faster. Refer to test_convergence.py to see the accuracy
# If it happens that user already have a pre-calculated value of Alm, she/he can provide it with the Alm_vals
# Skipping completely the computation.
def a_model_cpp(nu_nl, Dnu, a1, epsilon_nl, theta0, delta, ftype, l, Alm_vals=None):
	nu_nlm=[]
	if Alm_vals == None:
		el, em, Alm=Alm_cpp(l, theta0=theta0, delta=delta, ftype=ftype) # Array of m E [-l,l]
	else:
		Alm=Alm_vals[l-1] # Alm_vals contains a list of Alms for l=1,2,3 (index 0,1,2)
	#print('---------')
	for m in range(-l, l+1):	
		perturb_CF=nu_CF(nu_nl, Dnu, a1, l, m)
		perturb_AR=nu_nl*epsilon_nl*Alm[m+l]
		nu_nlm.append(nu_nl + perturb_CF + perturb_AR)
	acoefs=eval_acoefs(l, nu_nlm)
	return acoefs # returns all a-coeficients 

def a_model_interpol(nu_nl, Dnu, a1, epsilon_nl, theta0, delta0, l, interpolator_l1, interpolator_l2, interpolator_l3):
	nu_nlm=[]
	for m in range(-l, l+1):	
		if l==1:
			Alm=interpolator_l1[l+m](theta0, delta0)
		if l==2:
			Alm=interpolator_l2[l+m](theta0, delta0)
		if l==3:
			Alm=interpolator_l3[l+m](theta0, delta0)
		if l>=4:
			print("l>=4 is not implemented yet in a2_model_interpol")
			exit()
		if l<1:
			print("l<1 is not authorized in a2_model_interpol")
			exit()
		perturb_CF=nu_CF(nu_nl, Dnu, a1, l, m)
		perturb_AR=nu_nl*epsilon_nl*Alm
		nu_nlm.append(nu_nl + perturb_CF + perturb_AR)
	#print(nu_nlm)
	acoefs=eval_acoefs(l, nu_nlm)
	return acoefs

def nu_CF(nu_nl, Dnu, a1, l, m, a1_unit='nHz'):
	eta0=eta0_fct(Dnu=Dnu, rho=None)
	if a1_unit == 'nHz':
		return eta0*nu_nl * (a1*1e-9)**2 * Qlm(l,m)
	if a1_unit == 'microHz':
		return eta0*nu_nl * (a1*1e-6)**2 * Qlm(l,m)
	if a1_unit != 'nHz' and a1_unit != 'microHz':
		print('a1 must be provided either in nHz or in microHz')
		print('use the a1_unit argument of the nu_CF() function to set it properly')
		exit()

def nu_AR(nu_nl, epsilon_nl, theta0, delta, ftype, l):
	l,m,Alm=Alm_cpp(l, theta0=theta0, delta=delta, ftype=ftype)
	return nu_nl*epsilon_nl*Alm

def a2_CF(nu_nl, Dnu, a1, l):
	nu_nlm=[]
	for m in range(-l, l+1):
		perturb=nu_CF(nu_nl, Dnu, a1, l, m)
		nu_nlm.append(nu_nl + perturb)
	acoefs=eval_acoefs(l, nu_nlm)
	#print(nu_nlm)
	return acoefs[1] # returns only a2

def a2_AR(nu_nl, epsilon_nl, theta0, delta, ftype,l):
	nu_nlm=[]
	for m in range(-l, l+1):
		perturb=nu_AR(nu_nl, epsilon_nl, theta0, delta, ftype, l, m)
		nu_nlm.append(nu_nl + perturb)
	acoefs=eval_acoefs(l, nu_nlm)
	return acoefs[1] # returns only a2
