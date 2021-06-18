@read_ascii_kepler
pro ascii2sav, dir

	Nraws=1d6
	Ncols=3

	files=file_search(dir + '*.ascii')
	
	print, 'Looking for ascii files into the directory:'
	print, '     ' + dir
	for i=long(0), n_elements(files)-1 do begin
		b=byte(files[i])
		pslash=max(where(b eq 47)) ; symbol '/' (linux format)
		if pslash[0] eq -1 then begin
			pslash=max(where(b eq 47)) ; symbol '\' (windows format)
		endif
		pdot=max(where(b eq 46))  ; symbol '.'

		print, '['+strtrim(i,2) + '/' + strtrim(n_elements(files),2) +']' + ' Converting file: ' + strtrim(b[pslash+1:*], 2) + ' into ' + strtrim(b[0:pdot-1], 2) + '.sav'
		
		r=read_Ncolumns(files[i], Ncols, Nraws)
		freq=reform(r[0,*])
		spec_reg=reform(r[1,*])
		model=reform(r[2,*])
		
		
		save, freq, spec_reg, model, filename=strtrim(b[0:pdot-1], 2) + '.sav'
		;stop
	endfor

end
