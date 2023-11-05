'''
	Program that make a grid for all of the Alm terms
	Allows a faster computation in python later on
	The grid is required for the computation of the a-coefficients grid
	using show_aj_fct_theta_delta.py
'''
import numpy as np
import matplotlib.pyplot as plt
from activity import Alm_cpp
import time
from scipy import interpolate

def get_Alm(theta0, delta0, ftype, l, Alm_path='../cpp_prg/'):
	try:
		el, em, Alm=Alm_cpp(l, theta0=theta0, delta=delta0, ftype=ftype, Alm_path=Alm_path) # Array of m E [-l,l]
	except:
		print('ERROR ON ALM_CPP...')
		r=Alm_cpp(l, theta0=theta0, delta=delta0, ftype=ftype, Alm_path=Alm_path)
		print(r)
		exit()
	return Alm 

def make_Alm_grid(l=1, theta_min=0, theta_max=np.pi/2, delta_min=0, delta_max=np.pi/4, resol=np.pi/20, ftype='gauss', flatarray=False, Alm_path='../cpp_prg/'):

	Ntheta=int(np.ceil(theta_max-theta_min)/resol) # We maximise the number of point when rounding ==> resol is therefore the minimum required resolution
	Ndelta=int(np.ceil(delta_max-delta_min)/resol) # We maximise the number of point when rounding
	theta_axis=np.linspace(theta_min, theta_max, Ntheta)
	delta_axis=np.linspace(delta_min, delta_max, Ndelta)
	if flatarray == False:
		Alm=np.zeros((2*l+1, Ntheta,Ndelta))
		theta=theta_axis
		delta=delta_axis
	else:
		Alm=np.zeros((2*l+1, Ntheta*Ndelta))
		theta=np.zeros(Ntheta*Ndelta)
		delta=np.zeros(Ntheta*Ndelta)
	i=0 # Index in theta
	j=0 # Index in delta
	tot=0 # Linear flat inded for {Ntheta x Ndelta} space
	print('Number of data point on the theta axis:', Ntheta)
	print('Number of data point on the delta axis:', Ndelta)
	print('Resolution: ', resol)
	print('Theta range: ', '[', theta_min, theta_max, ']')
	print('Delta range: ', '[', delta_min, delta_max, ']')
	for theta0 in theta:
		print('theta0 =', theta0,  '   index:', str(i+1),   '/', str(Ntheta),'    (timestamp: ',str(time.time()),')')
		print('           index : [ 1 ,', Ndelta, ']')
		for delta0 in delta:
			#print('    (i, j)', '(',i, ',', j,')')
			#print('   tot: ', tot)
			v=get_Alm(theta0, delta0, ftype, l, Alm_path=Alm_path)
			#print('        ', v)
			if flatarray == False:
				Alm[:,i,j]=v
			else:
				theta[tot]=theta0
				delta[tot]=delta0
				Alm[:,tot]=v
			tot=tot+1
			j=j+1
		j=0
		i=i+1
	i=0

	resol_theta=theta[1]-theta[0] # Use the two first element of the regular grid to evaluate the effective resolution of the grid as it may differ (larger) than resol
	resol_delta=delta[1]-delta[0]
	return theta, delta, Alm, resol_theta, resol_delta

def do_Alm_grid(outfile='grid_Alm', l=1, resol=np.pi/20, ftype='gate', theta_min=0, theta_max=np.pi/2, delta_min=0, delta_max=np.pi/4, flatarray=False, Alm_path='../cpp_prg/'):
	m=np.linspace(-l, l, 2*l+1)
	print("             m      ")
	print("             ", m)
	theta,delta,Alm, resol_theta, resol_delta=make_Alm_grid(l=l, theta_min=theta_min, theta_max=theta_max, 
		delta_min=delta_min, delta_max=delta_max, 
		resol=resol, ftype=ftype, flatarray=flatarray, Alm_path=Alm_path)
	np.savez(outfile + "_" + str(l) + ".npz", theta=theta, delta=delta, Alm=Alm, l=l, m=m, resol_theta=resol_theta, resol_delta=resol_delta)
	grid=np.load(outfile + "_" + str(l) + ".npz")
	show_Alm_grids(grid, l, filerootjpg=outfile + '_', theta1=None, delta1=None)
	print('Grid Done')

def do_all_grids(dir_out=None, corefile='grid_Alm', resol=np.pi/20, ftype='gate', Alm_path='../cpp_prg/'):
	els=[1,2,3]
	for l in els:
		if dir_out != None:
			outfile=dir_out + '/' + corefile
		else:
			outfile=corefile
		if ftype == 'gate':
			delta_max=np.pi/4
		if ftype == 'triangle':
			delta_max=np.pi/2
		do_Alm_grid(outfile=outfile, l=l, resol=resol, ftype=ftype, theta_min=0, theta_max=np.pi/2, delta_min=0, delta_max=delta_max, flatarray=False, Alm_path=Alm_path)

def show_Alm_grids(grid, l, filerootjpg='grid_Alm_', theta1=None, delta1=None):
	#grid=np.load(gridfile)
	for m in range(-l, l+1):
		Alm_show=grid['Alm'][l+m,:,:]
		fig, ax = plt.subplots()
		ax.set_title("A"+str(l)+str(m))
		ax.set_xlabel('delta (rad)')
		ax.set_ylabel('theta (rad)')
		c = ax.pcolormesh(grid['delta'], grid['theta'],Alm_show, vmin=np.min(Alm_show), vmax=np.max(Alm_show), shading='auto')
		fig.colorbar(c, ax=ax)
#		if theta1 != None and delta1 != None:
#			ax.plot(theta1,delta1,marker='o',size=10)
		fig.savefig(filerootjpg + str(l) + "_" + str(m) + '.jpg')
		#plt.show()

def test_interpol_from_grid(gridfile, theta1, delta1, flatarray=False):

	if flatarray == True:
		print("Error: flatarray is not supported at the moment")
		exit()
	#
	grid=np.load(gridfile)
	resol_theta=grid['theta'][1] - grid['theta'][0]
	resol_delta=grid['delta'][1] - grid['delta'][0]
	theta_min=grid['theta'][0]
	delta_min=grid['delta'][0]
	Ntheta=len(grid['theta'])
	Ndelta=len(grid['delta'])
	#
	# Flatening of Alm
	l=grid['l']
	funcs=[]
	for m in range(-l,l+1):
		Alm_flat=[]
		for j in range(Ndelta):
			Alm_flat.append(grid['Alm'][l+m,:,j])
		funcs.append(interpolate.interp2d(grid['theta'], grid['delta'], Alm_flat, kind='cubic'))
	for m in range(-l,l+1):
		r=funcs[l+m](theta1, delta1)
		pos_theta=(theta1-theta_min)/resol_theta
		pos_delta=(delta1-delta_min)/resol_delta
		print(' m =', m)
		print('       Interpolated Alm   : ', r)
		print('       Nearest neightboors: ')
		print('       Fractional coords  : [pos(theta), pos(delta)] = [ ', pos_theta, ' , ', pos_delta, ' ]') 
		pos_theta_plus=int(np.ceil(pos_theta))
		pos_theta_minus=int(np.floor(pos_theta))
		pos_delta_plus=int(np.ceil(pos_delta))
		pos_delta_minus=int(np.floor(pos_delta))
		if pos_theta_plus == Ntheta:
			pos_theta_plus=-1
		if pos_delta_plus == Ndelta:
			pos_delta_plus=-1
		print('       Ceil and Floor  theta[',pos_theta_minus,']= ', grid['theta'][pos_theta_minus],  '      theta[',pos_theta_plus,']= ', grid['theta'][pos_theta_plus])
		print('                       delta[',pos_delta_minus,']= ', grid['delta'][pos_delta_minus],  '      delta[',pos_delta_plus,']= ', grid['delta'][pos_delta_plus])
		print('                       Alm for neightboors =')
		print('                               ',    grid['Alm'][l+m,pos_theta_minus,pos_delta_plus], '     ', grid['Alm'][l+m,pos_theta_plus,pos_delta_plus])
		print('                               ',    grid['Alm'][l+m,pos_theta_minus,pos_delta_minus], '     ', grid['Alm'][l+m,pos_theta_plus,pos_delta_minus])
		print('       ------      ')
	show_Alm_grids(grid, l, filerootjpg='gridplots_', theta1=theta1, delta1=delta1)
		#plt.show()