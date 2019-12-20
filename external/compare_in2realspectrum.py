import numpy
import scipy.io as scp
import matplotlib.pyplot as plt

# A simple function that read a .in file
def read_infile(infile):

	f=open(infile, 'r')

	alllines=f.readlines()
	tmp=alllines[0].strip()
	ID=float(tmp.split("=")[1])

	tmp=alllines[2].strip()
	tmp=tmp.split()
	Tobs=float(tmp[0])
	Cadence=float(tmp[1])

	# Modes section
	i0=3
	tmp=alllines[i0][0].strip()
	while tmp == "#":
		tmp=alllines[i0][0].strip()
		i0=i0+1
	iend=i0
	tmp=alllines[i0][0].strip()
	while tmp != "#":
		tmp=alllines[iend][0].strip()
		iend=iend+1

	i0=i0-1
	iend=iend-2
	Nrows=iend-i0

	Ncols=len(alllines[i0].split())
	mode_params=numpy.zeros((Nrows, Ncols))

	cpt=0
	for line in alllines[i0:iend]:
		tmp=line.split()
		for i in range(len(tmp)):
			#print("tmp[i]=", tmp[i])
			mode_params[cpt, i]=float(tmp[i])
		cpt=cpt+1
	f.close()

	# Noise section
	i0=iend+3
	
	Nrows=3
	Ncols=len(alllines[i0].split())
	noise_params=numpy.zeros((Nrows, Ncols))


	cpt=0
	for line in alllines[i0:i0+3]:
		tmp=line.split()
		#print("tmp = ", tmp)
		for i in range(len(tmp)):
			#print("tmp[i]=", tmp[i])
			noise_params[cpt, i]=float(tmp[i])
		cpt=cpt+1
	f.close()


	# Debug only
	'''
	print("ID=", ID)
	print("Tobs=", Tobs)
	print("Cadence=", Cadence)

	Nrows=len(mode_params[:, 0])
	Ncols=len(mode_params[0, :])
	print("Mode params:")
	print("		Nrows=", Nrows)
	print("		Ncols=", Ncols)
	for i in range(Nrows):
		stri=""
		for j in range(Ncols):
			stri=stri + "    " + str(mode_params[i, j])
		print(stri)

	Nrows=len(noise_params[:, 0])
	Ncols=len(noise_params[0, :])
	print("Noise params:")
	print("		Nrows=", Nrows)
	print("		Ncols=", Ncols)
	for i in range(Nrows):
		stri=""
		for j in range(Ncols):
			stri=stri + "    " + str(noise_params[i, j])
		print(stri)
	'''
	return mode_params, noise_params, ID, Tobs, Cadence

def main(infile, starfile, xmin=50, xmax=4000, starfile_format="sav"):

	mode_params, noise_params, ID, Tobs, Cadence=read_infile(infile)
	Nrows=len(mode_params[:,0])
	ind_el=0
	ind_freq=1
	ind_height=2
	ind_width=3

	amp=numpy.pi * mode_params[:, ind_height]*mode_params[:,ind_width]


	if starfile_format == "sav":
		data=scp.readsav(starfile)
		freq=data['freq']
		spec_reg=data['spec_reg']
	else:
		print("This function supports only IDL 'sav' files for the spectrum...")
		exit()


	plt.plot(freq, spec_reg)
	plt.xlim(xmin, xmax)
	color=[]
	marker=[]
	s=[]

	smax=100
	smin=1
	scaled_amp=(amp - min(amp))*(smax - smin)/(max(amp) - min(amp)) + smin # Used to create symbol of various size in function of the mode amplitude
#	print(scaled_amp)
#	exit()

	ymax=max(spec_reg[numpy.where(numpy.bitwise_and(freq <= xmax, freq>=xmin))])
	for i in range(Nrows):
		if mode_params[i, ind_el] == 0:
			color.append('black')
			#marker.append('o')
		if mode_params[i, ind_el] == 1:
			color.append('blue')
			#marker.append('o')
		if mode_params[i, ind_el] == 2:
			color.append('red')
			#marker.append('o')
		if mode_params[i, ind_el] == 3:
			color.append('brown')
			#marker.append('o')
		s.append(scaled_amp[i])
		plt.scatter(mode_params[i, ind_freq], ymax/4, s=s[i], c=color[i])#, marker=marker)
	plt.show()

infile='/Users/obenomar/Documents/GitHub/Spectra-Simulator-C/Data/Spectra_info/0000001.in'
starfile='/Users/obenomar/Documents/GitHub/Spectra-Simulator-C/test/9574283l.PE15/Input/9574283.sav'
xmin=200
xmax=700
main(infile, starfile, xmin=xmin, xmax=xmax)
#read_infile(infile)