# Evaluate the difference between outputs in python and outputs in C++
import numpy as np
from subprocess import Popen, PIPE
#import Alm as Alm_python
from activity import Alm as Alm_python
from activity import Alm_cpp

def test(ftype='gate'):
	Alm_path='../bin/'
	print("Test program for Alm computation")
	print("   This is to verify that the Alm computation in C++ matches the computation made in Python")
	print("   It is essential that this is the case because any analysis will relly on the exactness of the Alm term")
	print("   while for efficiency reason, it is prefereable to use C++ for doing this, it is unclear whether the integration")
	print("   that I use in C++ (Legendre integral in 2D) is accurate")

	#els=[1,2,3]
	els=[1]
	Ntests_theta=4
	Ntests_delta=10
	thetas=np.linspace(0, np.pi/2, Ntests_theta)
	deltas=np.linspace(0,np.pi/2, Ntests_delta)
	err=False
	err_limit=1e-3
	for l in els:
		for theta0 in thetas:
			for delta in deltas:
				l_c, m_c, Alm_c=Alm_cpp(l, theta0, delta, ftype, Alm_path=Alm_path)
				Alm_py=np.zeros(2*l+1)
				for m in range(-l, l+1):
					a=Alm_python(l,m, theta0=theta0, delta=delta, ftype=ftype)
					Alm_py[m+l]=a
				for i in range(len(m_c)):
					dif=100*(Alm_c[i]-Alm_py[i])/Alm_py[i]
					print("theta0 =", theta0,   "    delta =:", delta, "   l =", l_c[i], "   m =",m_c[i], "   Alm_cpp =", Alm_c[i], "   Alm_py =", Alm_py[i],   "    Dif (no unit)", Alm_c[i]-Alm_py[i], "    Delta (%) = ", dif)
					#if dif >= 0.1 and Alm_py[i] != 0:
					#	print("          WARNING: DIFFERENCE EXCEDING 1% WAS FOUND")
					#	err=True
					if Alm_c[i]-Alm_py[i] >= err_limit:
						print("          WARNING: DIFFERENCE EXCEDING {0:0.7f} WAS FOUND. This is {1:0.4f}% of the full scale [0,1]".format(err_limit, err_limit*100))
						err=True
	if err == True:
		print("Discrepancies between exceeding {} in absolute value between the C++ and the Python routine found. Please check their significance...".format(err_limit))

#test()