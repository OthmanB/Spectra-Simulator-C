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


; N: number of columns
; K: maximum number of lines
function read_Ncolumns, file,N, K, skip=skip, ref_N=ref_N, spectrum=spectrum

if n_elements(skip) eq 0 then skip=1 ; defaut we skip one line only
if n_elements(ref_N) eq 0 then ref_N=1 ; defaut we identify the zero non-used tab elements with column 1
if n_elements(spectrum) eq 0 then spectrum=1
openr, 3, file

	param=dblarr(K,N)
	a=''
	i=0d
      while EOF(3) ne 1 do begin
      	if i lt skip then readf,3,format='(q)'
        if i ge skip then begin
          	readf,3,a ; read data
          	uu=strsplit(a)
          	N_uu=N_elements(uu)-1
          	for j=0,N_uu-1 do begin
          		param(i,j)=float(strmid(a,uu(j),uu(j+1)-uu(j)-1))
          	endfor
			param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))
		endif
		i=i+1
      endwhile

close,3
param0=param
if ref_N ge 0 then begin
	test=where(param[*,ref_N] ne 0)
	param=param[test,*]
endif

print, 'END read'
if spectrum eq 0 then begin
	save, param, filename=file+'.sav'
endif else begin
	freq=dblarr(1, n_elements(param[*,0]))
	spec_reg=freq
	freq[0,*]=param[*,0]
	spec_reg[0,*]=param[*,1]
	save, freq, spec_reg, filename=file+'.sav'
endelse

return,param
end
