; Small program that look for ascii files in a given directory (dir)
; then make a plot of the spectra with proper smoothing and informations
@read_txt
pro idl_plots_spectra

	do_plots=1
	save_new=0 ; To save the data in a text file (useful because the first Harvey profile is not incorporated in my model)
	
	;dir_params='/home/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/Data_examples/M0.799/'
	;dir_spectra='/home/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/Data_examples/M0.799/Spectra_ascii/'

	;dir_params='/home/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/Data_examples/M0.915/'
	;dir_spectra='/home/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/Data_examples/M0.915/Spectra_ascii/'

    ; Linux:
	;dir_params='/home/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/Data_examples/M1.038/'
	;dir_spectra='/home/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/Data_examples/M1.038/Spectra_ascii/'
	;dir_out_ascii='/home/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/Data_examples/M1.038/data_rasha/'
	
	; Mac
	dir_params='/Users/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/Data_examples/M1.038/'
	dir_spectra='/Users/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/Data_examples/M1.038/Spectra_ascii/'
	;dir_out_ascii='/Users/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/Data_examples/M1.038/For_Rasha_samples/data_rasha/'
	
	;dir_params='/home/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/Data_examples/M1.096/'
	;dir_spectra='/home/obenomar/Dropbox/Temporary/Spectra-Simulator-Cpp_Wgrid/Data_examples/M1.096/Spectra_ascii/'

	
	dir_out=dir_params + 'Spectra_frames/'
	spawn, 'mkdir ' + dir_out
	spawn, 'mkdir ' + dir_out + 'eps/'
 
	file_params=file_search(dir_params + 'models_sample*.txt')
	if file_params eq '' then begin
		print, 'file with model parameters not found. '
		print, 'Searched syntax: ' + dir_params + 'models_sample*.txt'
		print, 'Check the filename/directory'
		print, 'The program will stop now'
		stop
	endif
	Ncols=61
	Nlines=100.
	skip=2
	d=read_txt(file_params,Ncols, Nlines, skip)
	
	Nlines_params=n_elements(d[0,*])
	col_logg=15
	col_numax=32
	col_dnu=31
	col_mass=1
	col_age=6

	Nlines_max=1000000.
	skip=1
	for i=0, Nlines_params-1 do begin
		fspec=dir_spectra + file_syntax(i, '', '.ascii')
		s=read_txt(fspec,2, Nlines_max, skip)

		H0=1000.
		tc=(1.d-6)/(5.d-6)
		p=2.
		Low_Harvey=H0/(1d + (s[0,*]*tc)^p)
		P1=randomu(seed, n_elements(Low_Harvey))
		P2=randomu(seed, n_elements(Low_Harvey))
		s[1,*]=s[1,*] + Low_Harvey*(P1^2 + P2^2)/2.

		resol=s[0,1]-s[0,0]
		scoef=d[col_dnu,i]*0.10/resol

		if do_plots eq 1 then begin
			set_plot,'ps'
			file_out=dir_out + 'frame' + string(i, f='(i05)')
			device, filename=file_out + '.eps', $
				/encap, /inch, /color;, xsize=10, ysize=10
			minx=1
			maxx=8000.
			miny=0.001
			maxy=1000.
			plot, s[0,*], s[1,*], /xst, /ylog, /xlog, ytitle='Power ' + textoidl('(ppm^2/\mu' + 'Hz)'), $
				xtitle='Frequency ' + textoidl('(\mu'+'Hz)'), charsize=1.75, $
				yr=[miny, maxy], xr=[minx, maxx], color=fsc_color('Black'), background=fsc_color('White'), /nodata, charthick=2, thick=2
			oplot, s[0,*], s[1,*], color=fsc_color('Grey')		
			oplot, s[0,*], smooth(s[1,*], scoef, /edge_truncate), color=fsc_color('Black')
			xyouts, 0.20, 0.30, 'Surface Gravity log(g)='+ string(d[col_logg,i], format='(f5.2)'),/normal, color=fsc_color('Black'), charsize=1.25, charthick=2
			xyouts, 0.20, 0.25, 'Star Age ='+ string(d[col_age,i], format='(f5.2)')  + ' (Gyrs)',/normal, color=fsc_color('Black'), charsize=1.25, charthick=2
			xyouts, 0.20, 0.20, 'Stellar Mass ='+ string(d[col_Mass,i], format='(f5.2)') + textoidl('M_{Sun}'),/normal, color=fsc_color('Black'), charsize=1.25, charthick=2
			;wait, 0.1
			device, /close

			spawn, 'convert -density 150x150 -flatten ' + $
				file_out  + '.eps '  + file_out + '.jpeg'

			spawn, 'mv ' + file_out + '.eps ' + dir_out + 'eps/' 
		endif
		if save_new eq 1 then begin
		    file_out_ascii= dir_out_ascii + string(i, f='(i05)') + '.ascii'
			openw, 3, file_out_ascii
				printf, 3, '#       freq (microHz)    spectrum (ppm2/microHz)'
				for j=long(0), n_elements(s[0,*])-1 do begin
					str=string(s[0,j], format='(f20.9)') + string(s[1,j], format='(f20.9)')
					printf, 3, str
				endfor
			close,3
		endif
	endfor

	;spawn, 'ffmpeg -r 20 -i ' + dir_out + 'frame%05d.jpeg ' + $
	;	' -qscale 1 ' + dir_out + 'movie.mp4'
	
end

function file_syntax, in, core, extension

	if in lt 10 then file_out=core+'000000'+strtrim(long(in),1)+ extension
	if in ge 10 AND in lt 100 then file_out=core+'00000'+strtrim(long(in),1)+ extension
	if in ge 100 AND in lt 1000 then file_out=core+'0000'+strtrim(long(in),1)+ extension
	if in ge 1000 AND in lt 10000 then file_out=core+'000'+strtrim(long(in),1)+ extension
	if in ge 10000 AND in lt 100000 then file_out=core+'00'+strtrim(long(in),1)+ extension
	if in ge 100000 AND in lt 1000000 then file_out=core+'0'+strtrim(long(in),1)+ extension
	if in ge 1000000 then file_out=core +strtrim(long(in),1)+'.sav'
	
return, file_out
end
