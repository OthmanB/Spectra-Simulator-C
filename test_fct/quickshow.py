import os
import numpy as np
import matplotlib.pyplot as plt
from zipfile import ZipFile

def read_datafile(file_in, ignore_syntax=False):
	'''
		Read a data file used for the TAMCMC analysis
		If the file is a zip file, then it is assummed that within it, there is a data file,
		that is going to be extracted first before reading
		file_in: The file to read
	'''
	extension=os.path.splitext(file_in)[1]
	if extension == '.zip':
		with ZipFile(file_in, 'r') as zip:
			#files=zip.printdir()
			files=zip.namelist()
			count=0
			for f in files:
				ext=os.path.splitext(f)[1]
				if ext == '.data':
					txt=zip.read(f).decode('utf-8')
					count=count+1
			if count > 1:
				print('Error: There was multiple data files within the zip files read by read_datafile()')
				print("       Debug required. File read:", file_in)
				exit()
	if extension == ".data":
		f=open(file_in, 'r')
		txt=f.read()
		f.close()
	if extension != ".data" and  extension != ".zip" and ignore_syntax==False:
		print("Error: The data file is neither a zip nor a data file")
		print("       Provided file: ", file_in)
		exit()
	if ignore_syntax == True and "txt" not in locals() and "txt" not in globals():
		f=open(file_in, 'r')
		txt=f.read()
		f.close()		
	data=txt.split('\n')
	Nlines=len(data)
	k=0
	header=[]
	for d in data:
		if d != "":
			key=d.strip()[0]
		else:
			key=""
		if key != "#" and key !='!' and key !='*' and key != "":
			line=d.split()
			if k == 0:
				output=np.zeros((Nlines, len(line)), dtype=float)
			output[k,:]=np.asarray(line, dtype=float)
			k=k+1
		else:
			header.append(d)
	if k < Nlines: # It means that there was a white line at the end of the file... need to remove it
		output=output[0:k,:]
	return output, header

def quickshow(x,y,m, xm=None, c=['red', 'orange', 'blue', 'cyan', 'purple'], do_loglog=False, fileout=None, xr=None, yr=None):
	try:
		nmodels=len(m[:,0])
	except:
		nmodels=1
		if xm is None:
			mnew=np.zeros((nmodels, len(x)))
			xm=x
		else:
			mnew=np.zeros((nmodels, len(xm)))
		mnew[0,:]=m
		m=mnew
	fig, ax= plt.subplots(2, 1, sharex=True, num=1, clear=True)
	if xr !=None:
		ax[0].set_xlim(xr)
		ax[1].set_xlim(xr)
	else:
		xr=[np.min(x), np.max(x)]

	if isinstance(yr, (list, tuple)) and len(yr) == 2:
		new_ylim = (yr[0], yr[1])  # Set the y-axis limits to the specified range
	elif yr == "max":
		# Calculate ymax and ymin based on the specified range
		pos = np.where(np.bitwise_and(x > xr[0], x < xr[1]))
		ymax = np.max(y[pos])
		ymin = np.min(y[pos])
		ax[0].set_ylim([ymin,ymax])
	elif yr ==None:	
		pass
	else:
		raise ValueError("yr should be a single value, a range [ymin, ymax], or 'max'")

	ax[0].plot(x, y, color='gray', label="Data")
	for j in range(nmodels):
		m_int=np.interp(x, xm, m[j, :])
		residuals = y / m_int
		ax[1].plot(x, residuals, color=c[j], label="Residuals " + str(j))
		if xm is None:
			ax[0].plot(x,m[j,:], color=c[j], label="Model " + str(j))
		else:
			ax[0].plot(xm,m[j,:], color=c[j], label="Model " + str(j))
	if do_loglog == True:
		ax[0].set_xscale('log')
		ax[0].set_yscale('log')
	ax[0].set_xlabel(r'Frequency $(\mu$Hz)')
	ax[0].set_ylabel(r'Power $(ppm/\mu$Hz)')
	ax[0].set_title('Quick Show Data vs Model')
	ax[0].legend()
	ax[1].plot(x, np.ones(len(x)), color='black', linestyle='--')
	ax[1].set_xlabel(r'Frequency $(\mu$Hz)')
	ax[1].set_ylabel('Residuals')
	plt.tight_layout()
	if fileout == None:
		plt.show()
	else:
		fig.savefig(fileout, dpi=300)

def do_show(data_file, fileout, do_loglog=False, xr=None, yr=None):
	'''
		data_file: data file full path
	'''
	current_dir=os.getcwd()
	print(" ... Reading data file...")
	data, hdr=read_datafile(data_file)
	print("...Showing the plot...")
	quickshow(data[:,0],data[:,1],data[:,2], xr=xr, yr=yr, c=['red', 'orange', 'blue', 'cyan', 'purple'], do_loglog=do_loglog, fileout=fileout)

def do_all_show(dir_in, dir_out, do_loglog=False):

	list_subdir = [subdir for subdir in os.listdir(dir_in) if not subdir.startswith('.')]
	# sort the list of subdir
	list_subdir.sort()
	# Loop over the list of subdir
	for process_name in list_subdir:
		print("    -----  Processing : ", process_name)
		filein=os.path.join(dir_in, process_name)
		fileout=os.path.join(dir_out, process_name+".jpg")
		do_show(filein, fileout, do_loglog=do_loglog)
		print("    -> File saved at : ", fileout)

print("  ----- Program to visualise outputs ----")
data_file = input("Enter the full path to the data file (3 column file): ")
fileout=None#"img.jpg"
do_loglog=True
do_show(data_file, fileout, do_loglog=False, xr=[1000,3000], yr="max")
