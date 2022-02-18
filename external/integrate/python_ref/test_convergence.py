# Evaluate the difference between outputs in python and outputs in C++ when using Alm
import numpy as np
from subprocess import Popen, PIPE
from activity import Alm as Alm_python

def Alm_cpp(l, theta0, delta, ftype, raw=False):
	process = Popen(["./Alm", str(l), str(theta0), str(delta), ftype], stdout=PIPE, stderr=PIPE)
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
	return -1, -1

def test():
	print("Test program for Alm computation")
	print("   This is to verify that the Alm computation in C++ matches the computation made in Python")
	print("   It is essential that this is the case because any analysis will relly on the exactness of the Alm term")
	print("   while for efficiency reason, it is prefereable to use C++ for doing this, it is unclear whether the integration")
	print("   that I use in C++ (Legendre integral in 2D) is accurate")

	ftype='gauss'
	els=[1,2,3]
	Ntests_theta=4
	Ntests_delta=20
	thetas=np.linspace(0, np.pi, Ntests_theta)
	deltas=np.linspace(0,np.pi/2, Ntests_delta)
	err=False
	for l in els:
		for theta0 in thetas:
			for delta in deltas:
				l_c, m_c, Alm_c=Alm_cpp(l, theta0, delta, ftype)
				Alm_py=np.zeros(2*l+1)
				for m in range(-l, l+1):
					a=Alm_python(l,m, theta0=theta0, delta=delta, ftype=ftype)
					Alm_py[m+l]=a
				for i in range(len(m_c)):
					dif=100*(Alm_c[i]-Alm_py[i])/Alm_py[i]
					print("theta0 =", theta0,   "    delta =:", delta, "   l =", l_c[i], "   m =",m_c[i], "   Alm_cpp =", Alm_c[i], "   Alm_py =", Alm_py[i],   "    Dif (no unit)", Alm_c[i]-Alm_py[i], "    Delta (%) = ", dif)
					if dif >= 0.1 and Alm_py[i] != 0:
						print("          WARNING: DIFFERENCE EXCEDING 1% WAS FOUND")
						err=True
	if err == True:
		print("Discrepancies between the C++ and the Python routine found. Please check their significance...")

#test()