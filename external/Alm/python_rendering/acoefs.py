import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from scipy.optimize import minimize
from scipy.special import legendre

def Hslm_Ritzoller1991(s,l, m):
	L=l*(l+1)
	if s == 0:
		Hsm=1
	if s == 1:
		Hsm=2*m
	if s == 2:
		Hsm=6*m**2 - 2*L
	if s == 3:
		Hsm=20*m**3 - 4*(3*L-1)*m
	if s == 4:
		Hsm=70*m**4 - 10*(6*L-5)*m**2 + 6*L*(L-2)
	if s == 5:
		Hsm=252*m**5 -140*(2*L-3)*m**3 + (20*L*(3*L-10) + 48)*m
	if s == 6:
		Hsm=924*m**6 - 420*m**4 * (3*L-7) + 84*m**2 *(5*L**2 - 25*L + 14) - 20*L*(L**2 - 8*L + 12)
	if s > 6:
		print("s>6 not supported")
		print("The program will return 0")
		return 0
	#beta_slm=(2*s+1)**(1./2) * Fs * Hsm
	return Hsm


def Pslm(s,l,m):
	# These take the Ritzoller1991 coefficients and normalise them by Pslm(l)=l
	# As per specified in Schou, JCD, Thompson 1994
	# so basically we solved Pslm(l)=l=c*Hslm and find c
	if s==0:
		Ps=l
	if s==1:
		Ps=m
	if s==2:
		Ps=(3*m**2 -l*(l+1))/(2*l-1)
	if s==3:
		Ps=(5*m**3 - (3*l*(l+1)-1)*m)/((l-1)*(2*l-1))
	if s==4:
		H=(35*m**4 - 5*(6*l*(l+1)-5)*m**2) + 3*l*(l+1)*(l*(l+1)-2)
		c=2*(l-1)*(2*l-1)*(2*l-3)
		Ps=H/c
	if s==5:
		H=Hslm_Ritzoller1991(s,l,m)
		c=8*(4*l**4 - 20*l**3 + 35*l**2 - 25*l + 6)
		Ps=H/c
	if s==6:
		H=Hslm_Ritzoller1991(s,l,m)
		c=64*l**5 - 480*l**4 + 1360*l**3 - 1800*l**2 + 1096*l - 240
		Ps=H/c
	if s>6:
		print("s>6 not supported")
		print("The program will return 0")
		return 0
	return Ps


def Pslm_Masao(s,l,m):
	'''
		Masao Used a formal programing language to derive the expressions of Pslm.
		The factorisation differed. This is to verify that the actual results are 
		the same as my hand-made normalisation
	'''
	L=l*(l+1)
	if s<=4:
		print("THERE IS NO DIFFERENCE WITH MASAO... IGNORED")
	else:
		if s==5:
			Ps=63*m**5 - 35*(2*L - 3)*m**3 + (15*L**2 - 50*L + 12)*m
			c=2*(2*l - 1)* (l - 1)*(2*l - 3) *(l - 2)
			Ps=Ps/c
		if s==6:
			Ps=231*m**6 - 105*(3*L-7)* m**4 + 21*(5*L**2 - 25*L + 14)*m**2 - 5*L*(L-6)*(L-2) 
			c=2*(2*l - 1)*(l - 1)*(2*l - 3)*(l - 2)*(2*l - 5)
			Ps=Ps/c
	return Ps

def Fls(l, s):
	Fs=scipy.special.factorial(2*l-s)/scipy.special.factorial(2*l+s+1)
	Fs=Fs**(1./2)
	return (2*s+1)**(1./2) * Fs

def Hm_s(smax, l, m):
	H0m=1
	H1m=2*m
	Hs_minus1_m=H0m
	Hsm=H1m
	for s in range(smax):
		Hs_plus1_m=2*(2*s+1)*m*Hsm - s*((2*l+1)**2 - s**2)*Hs_minus1_m/(s+1)
		Hs_minus1_m=Hsm
		Hsm=Hs_minus1_m
	return Hs_minus1_m

def eval_acoefs(l, nu_nls): # We expect nu_nls=[nu(-l), nu(-l+1), ... , nu[0], nu[1], ... nu(l)]
	# Function that gets the splitted frequencies of a given mode (n,l) and determines the analytical a-coefficients
	# from a1 to a6 and for l<=3
	if l ==0:
		print("There is no a1 coefficient for l=0")
	if l == 1: # ALL GOOD
		a1=(nu_nls[2]- nu_nls[0])/2
		a2=((nu_nls[0] + nu_nls[2])/2 - nu_nls[1])/3
		a3=0
		a4=0
		a5=0
		a6=0
	if l == 2: 
		Sn21=(nu_nls[1] + nu_nls[3])/2 - nu_nls[2]
		Sn22=(nu_nls[0] + nu_nls[4])/2 - nu_nls[2]
		Tn21=(nu_nls[3]-nu_nls[1])/2
		Tn22=(nu_nls[4]-nu_nls[0])/4
		# a1 Term:
		Num_a1=Tn21+4*Tn22 #T21-2*T22*Pslm(3,2,1)/Pslm(3,2,2)
		Den_a1=5 #Pslm(1,2,1) - Pslm(3,2,1)*Pslm(1,2,2)/Pslm(3,2,2)
		K_a1=1 #Pslm(1,2,2)/Pslm(3,2,2)
		a1=Num_a1/Den_a1
		#a1=(nu_nls[3] - nu_nls[1])/10 + (nu_nls[4]-nu_nls[0])/5
		# a3 term:
		a3=(Tn22 - Tn21)/5  # TO VERIFY
		# a2 term:
		a2=(2*Sn22 - Sn21)/7
		# a4 term:
		a4=(Sn22 - 4*Sn21)/70. 
		a5=0
		a6=0
	if l==3: # VERIFIED OK
		# Odds terms
		Tn31=(nu_nls[4]-nu_nls[2])/2
		Tn32=(nu_nls[5]-nu_nls[1])/4
		Tn33=(nu_nls[6]-nu_nls[0])/6
		# We have Tnlm = Sum [a_{2j-1} P_{2j-1}] with j=[1,M/2]
		# a1 term:
		a1=Tn31/14 + 2*Tn32/7 + 9*Tn33/14
		# a3 term: 
		a3=-Tn31/9 - 2*Tn32/9 + Tn33/3
		# a5 term: 
		a5 = Tn33/42 + 5*Tn31/126 - 4*Tn32/63
		# Even Terms
		Sn31=(nu_nls[2] + nu_nls[4])/2 - nu_nls[3]
		Sn32=(nu_nls[1] + nu_nls[5])/2 - nu_nls[3]
		Sn33=(nu_nls[0] + nu_nls[6])/2 - nu_nls[3]
		# We have Snlm = Sum[a2j (P2j - P2j(0))] with j=[1, M/2]
		'''
		Pi21=Pslm(2,3,1) - Pslm(2,3,0)		
		Pi22=Pslm(2,3,2) - Pslm(2,3,0)		
		Pi23=Pslm(2,3,3) - Pslm(2,3,0)		
		Pi41=Pslm(4,3,1) - Pslm(4,3,0)		
		Pi42=Pslm(4,3,2) - Pslm(4,3,0)		
		Pi43=Pslm(4,3,3) - Pslm(4,3,0)		
		Pi61=Pslm(6,3,1) - Pslm(6,3,0)		
		Pi62=Pslm(6,3,2) - Pslm(6,3,0)		
		Pi63=Pslm(6,3,3) - Pslm(6,3,0)		
		'''
		# a2 term:
		a2=(-15*Sn31 + 25*Sn33)/126
		# a4 term:
		a4=((Sn31 - 7*Sn32 + 3*Sn33)/91 )*13/11
		# a6 term:
		a6=(15*Sn31 - 6*Sn32 + Sn33)/1386
	if l>3:
		print('an for l>3 not implemented. Should you need it, better to solve this algorithmically using equation A3-6 from Schou, JCD, Thompson, 1994')
		print('The program will return 0 for a1,a2,a3,a4,a5,a6...')
		a1=0
		a2=0
		a3=0
		a4=0
		a5=0
		a6=0
	return a1,a2,a3,a4,a5,a6

def nunlm_from_acoefs(nunl0, l, a1=0, a2=0, a3=0,a4=0,a5=0,a6=0):
	nu_nlm=[]
	for m in range(-l, l+1):
		val=nunl0 + a1*Pslm(1,l,m) + a2*Pslm(2,l,m) 
		if l>=2:
			val=val + a3*Pslm(3,l,m) + a4*Pslm(4,l,m) 
		if l>=3:
			val=val + a5*Pslm(5,l,m) + a6*Pslm(6,l,m)
		nu_nlm.append(val)
	return nu_nlm

def Tnlm(nu_nlm, l):
	if l == 0:
		print(" Please enter a value of 0<l<4")
		return []
	if l == 1:
		Tn11=(nu_nlm[2]- nu_nlm[0])/2
		return [Tn11]
	if l == 2:
		Tn21=(nu_nlm[3]-nu_nlm[1])/2
		Tn22=(nu_nlm[4]-nu_nlm[0])/4
		return [Tn21, Tn22]
	if l == 3:
		Tn31=(nu_nlm[4]-nu_nlm[2])/2
		Tn32=(nu_nlm[5]-nu_nlm[1])/4
		Tn33=(nu_nlm[6]-nu_nlm[0])/6
		return [Tn31,Tn32, Tn33]
	if l > 3:
		print("Please enter a value of 0<l<4")
		return []

def Snlm(nu_nlm, l):
	if l == 0:
		print(" Please enter a value of 0<l<4")
		return []
	if l==1:   
		Sn11=(nu_nlm[0] + nu_nlm[2])/2 - nu_nlm[1]
		return [Sn11]
	if l==2:
		Sn21=(nu_nlm[1] + nu_nlm[3])/2 - nu_nlm[2]
		Sn22=(nu_nlm[0] + nu_nlm[4])/2 - nu_nlm[2]
		return [Sn21,Sn22]
	if l==3:
		Sn31=(nu_nlm[2] + nu_nlm[4])/2 - nu_nlm[3]
		Sn32=(nu_nlm[1] + nu_nlm[5])/2 - nu_nlm[3]
		Sn33=(nu_nlm[0] + nu_nlm[6])/2 - nu_nlm[3]
		return [Sn31,Sn32,Sn33]
	if l > 3:
		print("Please enter a value of 0<l<4")
		return []


def test_acoef(l=1, a1=1, a2=0, a3=0, a4=0, a5=0, a6=0, nu_nl=1000.):
	'''
		With this function you can see whether by putting some input
		a-coefficients, you do retrieve them after you measure directly
		the symetrical splittings (Tnlm) and asymetrical splitting (Snlm)
	'''
	nu_nlm=nunlm_from_acoefs(nu_nl, l, a1=a1, a2=a2, a3=a3,a4=a4,a5=a5,a6=a6)
	print("nu_nlm:", nu_nlm)

	print("  Tnlm : ")
	Tn=Tnlm(nu_nlm, l)
	for m in range(len(Tn)):
		print("Tn"+str(l)+str(m+1)+"  = ", Tn[m])		

	print(" Snlm : ")
	Sn=Snlm(nu_nlm, l)
	for m in range(len(Sn)):
		print("Sn"+str(l)+str(m+1)+"  = ", Sn[m])

	print("Input a-coefficients:")
	print("  anl =", [a1,a2,a3,a4,a5,a6])
	print("Retrieved a-coefficients from the function eval_acoefs():")
	anl=eval_acoefs(l, nu_nlm)
	print("  anl =", anl)