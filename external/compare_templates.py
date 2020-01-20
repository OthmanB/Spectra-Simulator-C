import numpy
import scipy.io as scp
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join

def read_templatefile(templatefile):
	try:
		f=open(templatefile, 'r')
		alllines=f.readlines()
		f.close()
	except:
		print("Error: Could not open the file ", templatefile)
		print("       Check that the file exists")
		print("       The program will exit now")
		exit()

	i0=0
	tmp=alllines[i0][0].strip()
	while tmp == "#": # Ignore comments in the header
		tmp=alllines[i0][0].strip()
		i0=i0+1

	i0=i0-1
	ID_ref=-1
	numax_ref=-1
	Dnu_ref=-1
	epsilon_ref=-1

	while tmp[0] != "#": # Look for keywords for numax, Dnu and ID
		tmp=alllines[i0].strip()
		keys=tmp.split("=")
		if keys[0] == "ID_ref":
			ID_ref=keys[1]
		if keys[0] == "numax_ref":
			numax_ref=float(keys[1])
		if keys[0] == "Dnu_ref":
			Dnu_ref=float(keys[1])
		if keys[0] == "epsilon_ref":
			epsilon_ref=float(keys[1])
		i0=i0+1

	if numax_ref == -1 or Dnu_ref == -1 or epsilon_ref == -1:
		print("Error: Could not find at least one of the keywords defining the global pulsation parameters")
		print("Check that the following keywords are present in the template file: ", templatefile)
		print("       The program will exit now")
		exit()
	if ID_ref == -1:
		print("Warning: ID_ref is not set")
		print("         This does not prevent the code to run, but may result in a more difficult tracking of the used reference mode profiles in the future")

	# Process the table
	Nrows=len(alllines[i0:])
	Ncols=len(alllines[i0].split())
	data=numpy.zeros((Nrows, Ncols))
	cpt=0
	for line in alllines[i0:]:
		tmp=line.split()
		#print("tmp = ", tmp)
		for i in range(len(tmp)):
			#print("tmp[i]=", tmp[i])
			data[cpt, i]=float(tmp[i])
		cpt=cpt+1

	return ID_ref, numax_ref, Dnu_ref, epsilon_ref, data

def show_all(templatedir):

	templatefiles=[f for f in listdir(templatedir) if isfile(join(templatedir, f))]

	ymax=-1
	for file in templatefiles:
		ID_ref, numax_ref, Dnu_ref, epsilon_ref, data=read_templatefile(file)
		freq=data[:,0]
		height=data[:,1]
		plt.plot(freq, height)
		if ymax < max(height):
			ymax=max(height)
	
	plt.show()

	ymax=-1
	for file in templatefiles:
		ID_ref, numax_ref, Dnu_ref, epsilon_ref, data=read_templatefile(file)
		freq=data[:,0]
		width=data[:,2]
		plt.plot(freq, width)
		if ymax < max(width):
			ymax=max(width)
	
	plt.show()


def compare_two(templatedir, file1, file2):
	ymax=-1
	ID_ref1, numax_ref1, Dnu_ref1, epsilon_ref1, data1=read_templatefile(templatedir+file1)
	freq1=data1[:,0]
	height1=data1[:,1]
	width1=data1[:2]

	ID_ref2, numax_ref2, Dnu_ref2, epsilon_ref2, data2=read_templatefile(templatedir+file2)
	freq2=data2[:,0]
	height2=data2[:,1]
	width2=data2[:2]
	
	plt.plot(freq1, height1)
	plt.plot(freq2, height2)
	plt.show()

templatedir='/Users/obenomar/Documents/GitHub/Spectra-Simulator-C/Configurations/templates/'
#show_all(templatedir)


file_ref='Sun.template'
file2='Sun-mod.template'
compare_two(templatedir, file_ref, file2)