import numpy
import scipy.io as scp
import matplotlib.pyplot as plt

def smooth(x,y,window_len=11,window='flat'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    #if x.ndim != 1:
    #    raise ValueError, "smooth only accepts 1 dimension arrays."
    #
    #if x.size < window_len:
    #    raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return y


    #if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
    #    raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[y[window_len-1:0:-1],y,y[-2:-window_len-1:-1]]
    
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    r=numpy.convolve(w/w.sum(),s,mode='valid')
    Nr=len(r)
    #x=x[int(window_len/2):Nx-int(window_len/2)-1]
    rout=r[int(window_len/2):Nr-int(window_len/2)]
    if len(rout) != len(x):
    	rout=r[int(window_len/2):Nr-int(window_len/2)-1]
    if len(rout) != len(x):
    	print("error in smoothing. Debug required")
    	exit()
    return x,rout



def read_datfile(asciifile):
	f=open(asciifile, 'r')
	alllines=f.readlines()
	f.close()

	i0=0
	tmp=alllines[i0][0].strip()
	while tmp == "#":
		tmp=alllines[i0][0].strip()
		i0=i0+1

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

	return data

# A simple function that read a .in file
def read_infile(infile):

	f=open(infile, 'r')

	alllines=f.readlines()
	f.close()

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

def main(infile, starfile, asciifile=None, xmin=50, xmax=4000, starfile_format="sav"):

	mode_params, noise_params, ID, Tobs, Cadence=read_infile(infile)
	Nrows=len(mode_params[:,0])
	ind_el=0
	ind_freq=1
	ind_height=2
	ind_width=3

	if asciifile != None:
		data_simu=read_datfile(asciifile)

	amp=numpy.pi * mode_params[:, ind_height]*mode_params[:,ind_width]


	if starfile_format == "sav":
		data=scp.readsav(starfile)
		freq=data['freq']
		spec_reg=data['spec_reg']
	else:
		print("This function supports only IDL 'sav' files for the spectrum...")
		exit()

	#spec_reg
	f,s=smooth(freq, spec_reg,window_len=11,window='flat')
	plt.plot(f, s)
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
		if asciifile != None:
			plt.plot(data_simu[:,0], data_simu[:,2])#+ 1.5)
	plt.show()

#infile='/Users/obenomar/Documents/GitHub/Spectra-Simulator-C/Data/Spectra_info/0000001.in'

#infile='/Users/obenomar/Documents/GitHub/Spectra-Simulator-C/test/9574283l.Simu/Spectra_info/0000001.in'
#datafile='/Users/obenomar/Documents/GitHub/Spectra-Simulator-C/test/9574283l.Simu/Spectra_ascii/0000001.ascii'
#starfile='/Users/obenomar/Documents/GitHub/Spectra-Simulator-C/test/9574283l.PE15/Input/9574283.sav'

infile='/Users/obenomar/Documents/GitHub/Spectra-Simulator-C/test/8751420.Simu/Spectra_info/0000001.in'
datafile='/Users/obenomar/Documents/GitHub/Spectra-Simulator-C/test/8751420.Simu/Spectra_ascii/0000001.ascii'
starfile='/Users/obenomar/Documents/GitHub/Spectra-Simulator-C/test/8751420.PE15/Input/8751420.sav'

xmin=200
xmax=700
#main(infile, starfile, xmin=xmin, xmax=xmax)
main(infile, starfile, asciifile=datafile, xmin=xmin, xmax=xmax)
#read_infile(infile)