'''
	Module that allows you to show a butterfly diagram for the Sun
	It can also show superimposed, models of the Alm filter
'''
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from datetime import date
from sunpy.coordinates.sun import carrington_rotation_time
from sunpy.coordinates.sun import carrington_rotation_number
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy.signal import medfilt
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from activity import gate_filter, gauss_filter, gauss_filter_cte, triangle_filter 

def read_dailyarea(file, valid_only=True, date_format='original'):
	'''
		Read the daily average number of spots from 1874 (until 2016 as the date this code is writen - 02/06/2022)
		as per provided at https://solarscience.msfc.nasa.gov/greenwch.shtml
		file: input data file
		valid_only: If True, remove negative (-1) counts corresponding to days without observations
		date_format: Set the output format of the data
			- 'original': Keep the (YY, MM, DD) data for the date as it is
			- 'date': convert the 3 columns giving the time (YY-MM-DD) into a single input of type datetime. 
					  Note that the output here is a list for date and a numpy array for the other data (Total, North, South)
		    - 'timestamp': convert the 3 columns giving the time (YY-MM-DD) into a single timestamp as per defined
		    			   in the datetime python package
	'''
	# Read the file
	f=open(file, 'r')
	txt=f.read()
	f.close()
	txt=txt.split('\n')
	# Separate the first line (header) from the raw data
	header=txt[0]
	raw_data=txt[1:]
	# Parse the raw data into a numpy table
	Ncols=len(raw_data[0].split())
	Nlines=len(raw_data)
	data=np.zeros((Nlines, Ncols))
	cpt=0
	for i in range(Nlines):
		line=raw_data[i].split()
		if line != '' and line !=[]:
			if valid_only == False:
				data[i, :]=line
			else:
				data[cpt, :]=line 
				cpt=cpt+1
	if valid_only == True:
		data=data[0:cpt,:]
		Nlines=cpt-1
	# Handling convert_to_date
	if date_format == 'original':
		return data
	if date_format == 'date':
		date=[]
		outputs=np.zeros((Nlines, Ncols-3))
		for i in range(Nlines):
			date.append(datetime(int(data[i,0]), int(data[i,1]), int(data[i,2])))
			outputs[i,:]=data[i,3:]
		return date, outputs
	if date_format == 'timestamp':
		outputs=np.zeros((Nlines, Ncols-2))
		for i in range(Nlines):
			outputs[i,0]=datetime.timestamp(datetime(int(data[i,0]), int(data[i,1]), int(data[i,2])))
			outputs[i,1:]=data[i,3:]
		return outputs

def read_butterflydata(file, date_format='carrington'):
	# Read the file
	f=open(file, 'r')
	txt=f.read()
	f.close()
	txt=txt.split('\n')
	# Separate the first line (header) from the raw data
	header=txt[0]
	raw_data=txt[1:]
	Ncols=len(raw_data[0].split())
	Nlines=len(raw_data)
	Blocksize=6 # A data block is 1 Carington number, followed by a data block of 5 lines
	Nblocks=int(Nlines/Blocksize) # Determine the Number of blocks
	# Read the Carrington numbers only first
	carrington=np.zeros(Nblocks, dtype=float)
	cpt=0
	for i in range(0,Nlines, Blocksize):
		carrington[cpt]=raw_data[i]
		cpt=cpt+1
	# Read the data within a block. The number of records within a block is 50,
	# distributed within 5 lines of 10 entry each. 
	# Those 50 entries give the number of spots as a function sin(latitude)
	spot_area=np.zeros((Nblocks,50))
	cpt=0
	for i in range(1, Nlines, Blocksize): # For each block
		for j in range(Blocksize-1): # Read each block lines
			line=raw_data[i+j].split(',')[0:-1] # We remove the empty '' at the end due to the last ','
			if line != '' and line !=[]:
				spot_area[cpt, j*len(line):(j+1)*len(line)]=line
		cpt=cpt+1
	# Compute the latitude array in degree: Each 50 records from the file is a sin(latitude) 
	latitude=np.arcsin(np.linspace(-1, 1, 50))*180./np.pi

	# Give the date as per specified by the user
	if date_format == 'carrington':
		time=carrington
		passed=True
	if date_format == 'date': 
		# Convert the Carrington number to an approximate date in datetime package format
		# Exact time is ignored so that we might have an error of 12h
		time=[]
		for c in carrington:
			t=date.fromisoformat(carrington_rotation_time(c).value.split()[0])
			print('     original: ', c, '    value: ', carrington_rotation_time(c).value, '    time:',t)
			time.append(t)
		passed=True
	if date_format == 'timestamp':
		# Convert the Carrington number to an exact timestamp as per defined in the datetime package
		time=[]
		for c in carrington:
			t=datetime.fromisoformat(carrington_rotation_time(c).value).timestamp()
			print('     original: ', c, '    value: ', carrington_rotation_time(c).value, '    time:',t)
			time.append(t)
		time=np.asarray(time)
		passed=True
	if passed != True:
		print("Error: Invalid provided date_format. Set it to either 'carrington', 'date' or 'timestamp'")
	return time, latitude, spot_area.transpose()

def select_dates(date_interval, carrington_time, spot_area):
	'''
		Function that retrieve a provided date interval
		date_interval: Array of string with a date in the ISO format ('YY-MM-DD')
		carrington_time: Time as read in the Greenwich table for the butterfly diagram: In Carrington Units
		spot_area: A 2D numpy array as per returned by read_butterflydata()
	'''
	# Convert the date interval into a Carrington time
	tmin=carrington_rotation_number(datetime.fromisoformat(date_interval[0]))
	tmax=carrington_rotation_number(datetime.fromisoformat(date_interval[1]))
	posOK=np.where(np.bitwise_and(carrington_time >= tmin, carrington_time <= tmax))
	return carrington_time[posOK[0]], spot_area[:,posOK[0]]

def select_latitudes(latitude_interval, latitude, spot_area):
	'''
		Function that retrieve a provided latitude interval
		latitude_interval: Array with two values: A minimum and maximum latitude
		latitude: The full latitude range associated to spot_area
		spot_area: A 2D numpy array as per returned by read_butterflydata()
	'''
	posOK=np.where(np.bitwise_and(latitude >= latitude_interval[0], latitude <= latitude_interval[1]))
	return latitude[posOK[0]], spot_area[posOK[0],:]


def get_year_frac(date_in):
	try:
		start = date(date_in.year, 1, 1).toordinal()
	except: # In case the date was actually provided as a string an not as a datetime object... try to convert it and retry
		date_in=datetime.fromisoformat(date_in)
		start = date(date_in.year, 1, 1).toordinal()
	year_length = date(date_in.year+1, 1, 1).toordinal() - start
	return date_in.year + float(date_in.toordinal() - start) / year_length

def carrington2fracyear(carrington_time):
	d=datetime.fromisoformat(carrington_rotation_time(int(carrington_time)).value)
	return get_year_frac(d)

def do_diagram(ax, date_interval, latitude, spot_area, latitude_range=[-90, 90], timezones = [], text_index=None):
	unit=1e-3 
	vmin=0
	vmax=np.max(spot_area*unit)/10
	# Convert map in % of the Hemispheric area. Raw units are in 10^-6 of Hemisphere, so we need 
	# to multiply them by 10^-3 to get % of Hemisphere
	l, s = select_latitudes(latitude_range, latitude, spot_area*unit)
	im = ax.imshow(s, aspect='auto', cmap='Reds', vmin=vmin, vmax=vmax, extent = [date_interval[0] , date_interval[1], latitude_range[0] , latitude_range[1]])
	ax.set_xlabel('Year', fontsize=14)
	ax.set_ylabel('Latitude (deg)', fontsize=14)
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)
	ax.axhline(0, linestyle='--', color='black')
	cbaxes = inset_axes(ax, width="20%", height="2%", loc='upper right') 
	bar1 = plt.colorbar(im, orientation='horizontal', cax=cbaxes, ticks=[round(vmin,2), round(vmax,2)])
	bar1.ax.tick_params(labelsize=9)
	if text_index != None:
		ax.annotate(text_index, xy=(0.88, 0.05), xycoords=ax.transAxes, fontsize=18)
	# If colored zones are requested, handle them
	# One expects zones to be a list of lists. The lower level list
	# must be of the format [min_date, max_date, color_name]
	# min and max dates are in string of the ISO format 'YY-MM-DD' or datetime objects
	if timezones != []:
		alpha=0.4
		for z in timezones:
			rectangle = Rectangle((get_year_frac(z[0]), latitude_range[0]), get_year_frac(z[1])-get_year_frac(z[0]), latitude_range[1]- latitude_range[0], edgecolor=None, facecolor=z[2], linewidth=None, alpha=alpha)
			ax.add_patch(rectangle)
 
def do_projection(ax, date_filters, carrington_time, latitude, spot_area, rotate=False, latitude_range=[-45,45], smooth_lvl=1, norm=None, text_index=None):
	Nzones=len(date_filters)
	Nlatitudes=len(spot_area[:,0])
	norms=[] # To return the y-axis maximum heights of the collapsogram
	for i in range(Nzones):
		if len(date_filters[i]) == 2:
			dates=date_filters[i]
		else:
			dates=date_filters[i][0:2]
		c, s=select_dates(dates, carrington_time, spot_area)
		collapsogram=np.zeros(Nlatitudes)
		# by default the collapsogram uses the 'mean' for the normalisation 
		# and convert into % of the Hemispheric area raw units are in 10^-6 of Hemisphere, so we need 
		# to multiply them by 10^3 to get % of Hemisphere
		if norm == None:
			norm=len(s[0,:])*1e3
		#	
		for j in range(Nlatitudes):
			collapsogram[j]=np.sum(s[j,:])/norm # Collapse the time axis to have the total spot area over period
		# Add indication of the median
		pos_north=np.where(latitude > 0)
		pos_south=np.where(latitude < 0)
		med_north=np.average(latitude[pos_north], weights=collapsogram[pos_north])
		med_south=np.average(latitude[pos_south], weights=collapsogram[pos_south])
		collapsed = medfilt(collapsogram,smooth_lvl)
		norms.append(np.max(collapsed))
		if rotate == False:
			ax.plot(latitude, collapsed, label=date_filters[i][0] + ' to ' + date_filters[i][1])
			ax.set_xlim(latitude_range)
			ax.set_xlabel('Latitude (deg)', fontsize=14)
			ax.set_ylabel('Spot Area '+ r'($\%$ Hemisphere)', fontsize=14)	
			ax.axvline(0, linestyle='--', color='black')	
			if date_filters[i][3] != False:
				ax.axvline(med_north, linestyle='-', color=date_filters[i][2], xmin=0.9, xmax=1)
				ax.axvline(med_south, linestyle='-', color=date_filters[i][2], xmin=0.9, xmax=1)
		else:
			ax.plot(collapsed, -latitude, color=date_filters[i][2],label=date_filters[i][0] + ' to ' + date_filters[i][1])
			ax.set_ylim(latitude_range)	
			ax.set_xlabel('Area '+ r'($\%$)', fontsize=14)
			ax.tick_params(axis='x', labelsize=14)
			ax.axhline(0, linestyle='--', color='black')
			if date_filters[i][3] != False:
				ax.axhline(med_north, linestyle='-', color=date_filters[i][2], xmin=0.9, xmax=1, linewidth=2)
				ax.axhline(med_south, linestyle='-', color=date_filters[i][2], xmin=0.9, xmax=1, linewidth=2)
				ax.annotate("{0:.0f}".format(med_north), xy=(0.188, med_north), fontsize=8, va="center")
				if i ==0:
					ax.annotate("{0:.0f}".format(med_south), xy=(0.188	, med_south), fontsize=8, va="center")
			ax.get_yaxis().set_visible(False)
	# Handling legends
	ax.legend(fontsize=5, loc='upper left')	
	if text_index != None:
		ax.annotate(text_index, xy=(0.65, 0.05), xycoords=ax.transAxes, fontsize=18)
	#
	return norms

#import numpy as np
#import matplotlib.pyplot as plt
#fig, ax= plt.subplots(1, figsize=(12, 12))
#model_params=['triangle', 75, 10]
def show_spot_model(ax, model_params, norm=1, rotate=True, color='Black', linestyle='-'):
	'''
		The core function that plots the model of the spots
		ax: The plot zone
		model_params: a 3-element lists with:
			[0] the type of model (triange, gate, gauss)
			[1] the theta0 parameter of the model (COLATITUDE)
			[2] the delta parameter of the model (ACTIVE ZONE EXTENSION)
		latitude_range: The range of latitudes to be computed and showed
		norm: The maximum height of the model. Used to scale when superimposing
		on another plot
	'''
	cte=180./np.pi # To convert the native pi units into degrees
	theta_min=0
	theta_max=np.pi
	theta=np.linspace(theta_min, theta_max, 1000)
	if model_params[0] == 'triangle':
		F=triangle_filter(theta, model_params[1]/cte, model_params[2]/cte)
	if model_params[0] == 'gate':
		F=gate_filter(theta, model_params[1]/cte, model_params[2]/cte)
	if model_params[0] == 'gauss':
		F=gauss_filter(theta, model_params[1]/cte, model_params[2]/cte)
		F=F/gauss_filter_cte(model_params[1]/cte, model_params[2]/cte)
	theta_colat=90 - theta*cte # Convert the colatitudes in latitudes
	if rotate == False:
		ax.plot(theta_colat, F*norm, color=color, linestyle=linestyle)
	else:
		ax.plot(F*norm, -theta_colat, color=color, linestyle=linestyle)
	#print(theta_colat)
	#ax.axvline(x=(90 - model_params[1])/cte, color='green', linestyle='-.')

def show_diagram(filein=None, fileout=None, add_model=['triangle', 75, 10]):
	'''
		Generate an image of the butterfly diagram by reading the data stored in the filein
		Optionally, a spot-area model can be overplotted
		filein: The input file for the Butterfly diagram. If not provided, look in the relative path expected 
				from the paper Benomar+2022 (2023).
		fileout: The output image file. If not provided, look in the relative path expected 
				from the paper Benomar+2022 (2023).
		add_model: Add a model for the shape of the activity zone in the projected butterfly diagram
				   Type of profile available: 'gate' (a gate function), 'gauss' (a gaussian function), 'triangle'
				   If set to an empty list [], this function is deactivated.
	'''
	cwd = os.getcwd()
	if filein == None:
		filein = cwd + '/../../Data/External_data/Greenwich-data/butterflydata.txt'
	if fileout == None:
		fileout = cwd + '/../../Data/Figures_publish/Fig2-butterfly.jpg'
	#filein='/Users/obenomar/Work/tmp/test_a2AR/tmp/External-data/Greenwich-data/butterflydata.txt' #-1985-2016.txt'
	#
	#fileout='/Users/obenomar/Work/tmp/test_a2AR/tmp/External-data/Greenwich-data/butterfly-2006.5-2010.5.jpg'

	date_filter_butterfly=['1985-01-01', '2022-01-01']
	#
	date_filters_projection=[]
	date_filters_projection.append(['1986-01-01', '2016-01-01', 'blue', True])
	#date_filters_projection.append(['1988-01-01', '1992-01-01', 'purple', True])
	#date_filters_projection.append(['1999-01-01', '2002-01-01', 'blue',  True])
	#date_filters_projection.append(['2006-01-01','2009-01-01',   'orange', False])
	#
	carrington_time, latitude, spot_area=read_butterflydata(filein, date_format='carrington')
	carrington_time, spot_area=select_dates(date_filter_butterfly, carrington_time, spot_area)

	# Compute the min and max date of the data set, assuming it was provided in Carrington units
	date_interval=[carrington2fracyear(carrington_time[0]), carrington2fracyear(carrington_time[-1])]
	print('Initial Carrington_time =', int(carrington_time[0]),  '    Date: ', date_interval[0])
	print('Final   Carrington_time =',int(carrington_time[-1]), '     Date: ', date_interval[1])
	fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]})
	#fig, ax= plt.subplots(1, figsize=(12, 12))
	do_diagram(ax[0], date_interval, latitude, spot_area, latitude_range=[-60, 60], timezones=date_filters_projection, text_index='(a)')
	norm_collapsogram=do_projection(ax[1], date_filters_projection, carrington_time, latitude, spot_area, rotate=True, latitude_range=[-60, 60], text_index='(b)')
	if add_model != []:
		show_spot_model(ax[1], add_model, norm=np.max(norm_collapsogram), rotate=True, color='Black', linestyle='--')
	fig.tight_layout()
	plt.savefig(fileout, dpi=300)

	print('file save at: ', fileout)
	
filein='../data/Greenwich-data/butterflydata.txt'
fileout='../data/outputs/butterfly.jpg'
show_diagram(filein=filein, fileout=fileout, add_model=['triangle', 75, 30])