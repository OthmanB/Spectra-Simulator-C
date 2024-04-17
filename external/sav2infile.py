#import numpy as np
import scipy.io as scp

def sav2in_global(file_in, file_out='out.in'):
	params=scp.readsav(file_in)

	lmax=int(params['parameters_length'][1]) + 1
	Nmax=len(params['stat_synthese_freq'][:,0,0])

	if len(params['parameters_length']) == 9:
		fitversion='idl'
	else:
		fitversion='cpp'

	b=0
	alfa=0
	if fitversion == "idl":
		pos_a1=int(sum(params['parameters_length'][0:3]))
		pos_eta=int(sum(params['parameters_length'][0:3]) +1)
		pos_a3=int(sum(params['parameters_length'][0:3]) + 2)
		pos_asym=int(sum(params['parameters_length'][0:3]) + 5)
		pos_inc=int(sum(params['parameters_length']) -1)

		pos_h1=int(sum(params['parameters_length'][0:5]))
		Nharvey=int( (params['parameters_length'][5]-1)/3)
	else:
		pos_a1=int(sum(params['parameters_length'][0:6]))
		pos_eta=int(sum(params['parameters_length'][0:6]) +1)
		pos_a3=int(sum(params['parameters_length'][0:6]) + 2)
		pos_asym=int(sum(params['parameters_length'][0:6]) + 5)
		pos_inc=int(sum(params['parameters_length']) -1)

		pos_h1=int(sum(params['parameters_length'][0:8]))
		Nharvey=int( (params['parameters_length'][8]-1)/3)


	a1=params['stat_synthese_unsorted'][pos_a1][3]
	eta=params['stat_synthese_unsorted'][pos_eta][3]
	a3=params['stat_synthese_unsorted'][pos_a3][3]
	asym=params['stat_synthese_unsorted'][pos_asym][3]
	inc=params['stat_synthese_unsorted'][pos_inc][3]

	fout=open(file_out, 'w')

	header="# This file contains the results converted from an IDL sav file for the star with unknown ID: \n"
	fout.write(header)
	ID="-1 \n"
	fout.write(ID)

	header="# Spectrum parameters. Observation duration (days) / Cadence (seconds) \n"
	fout.write(header) 
	header=" -1    -1 \n"
	fout.write(header) 
	header="# Input mode parameters. degree / freq / H / W / splitting a1 / eta / a3 / b  / alfa  / beta asym /inclination \n"
	fout.write(header)
	
	for el in range(0, lmax):
		for en in range(0, Nmax):
			nu=params['stat_synthese_freq'][en, el, 3]
			h=params['stat_synthese_height'][en, el, 3]
			try:
				w=params['stat_synthese_width'][en, el, 3]
			except:
				w=params['stat_synthese_largeur'][en, el, 3]
			if (nu > 0) and (h>0) and (w>0):
				stri='{0}   {1:.5f}   {2:.5f}   {3:.5f}   {4:.4f}   {5:.10f}    {6:.5f}    {7}    {8}   {9:.5f}    {10:.3f}'.format(el, nu, h, w, a1, eta, a3, b, alfa, asym, inc)
				fout.write(stri)
				fout.write('\n')

	header="# Input noise parameters. H0, tau_0, p0 / H1, tau_1, p1 / N0. Set at -1 if not used. \n"
	fout.write(header)
	cpt=0
	for n in range(Nharvey):
		h_n=params['stat_synthese_unsorted'][pos_h1+cpt][3]
		tau_n=params['stat_synthese_unsorted'][pos_h1+cpt+1][3]
		p_n=params['stat_synthese_unsorted'][pos_h1+cpt+2][3]
		if h_n <= 0 or tau_n <=0 or p_n <=0:
			h_n=-1
			tau_n=-1
			p_n=-1

		stri='{0:.5f}   {1:.5f}   {2:.5f}'.format(h_n, tau_n, p_n)
		stri=stri + "\n"
		fout.write(stri)
		cpt=cpt+3
	pos_N0=int(pos_h1+3*Nharvey)
	N0=params['stat_synthese_unsorted'][pos_N0][3]
	fout.write('{0:.5f}'.format(N0))
	fout.close()
	return params


def sav2in(file_in, file_out, fittype='global'):
	if fittype=='global':
		sav2in_global(file_in, file_out)
	else:
		print("Warning: Conversion from sav file is only available for fittype='global'")
		print("         Cannot proceed further. The program will exit now")
		exit()

def main(file_in, file_out, fittype):
	sav2in(file_in, file_out, fittype=fittype)
	print('sav file converted successfully')
	exit()