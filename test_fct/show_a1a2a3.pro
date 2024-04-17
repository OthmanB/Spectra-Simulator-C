@fsc_color
pro show_allinfiles, dir_in

	dir_in='/Users/obenomar/tmp/Simulator/Spectra-Simulator-C-1.0.0/test_fct/cfg_plots/build_l_mode_a1a2a3/'
	files=file_search(dir_in + '*.in')
	for i=0, n_elements(files)-1 do begin
		print, '[' + strtrim(i,2) + '] Processing: ' + files[i] + '...'
		show_a1a2a3,files[i]
	endfor
end
;Visualisation for the cfg of test_lorentzian_a1a2a3
pro show_a1a2a3, cfg_file, c_prg_path=c_prg_path
	model_name='build_l_mode_a1a2a3'
	;cfg_file='model_l1-1.in'
	;cfg_file='model_l1-2.in'
	;cfg_file='model_l1-3.in'
	;cfg_file='model_l1-4.in'
	;cfg_file='model_l2-1.in'
	;cfg_file='model_l2-2.in'
	;cfg_file='model_l2-3.in'
	;cfg_file='model_l2-4.in'

	if n_elements(c_prg_path) eq '' then begin
		c_prg_path='test_lorentzians'
	endif
	f=file_search(cfg_file)
	if f[0] ne '' then begin
		inputs=read_cfg_file(cfg_file)
	endif else begin
		print, 'Could not found the specified configuration file : ', cfg_file
		print, 'Check the location of the file'
		print, 'The program will stop now'
		stop		
	endelse

	f=file_search(c_prg_path)
	if f[0] ne '' then begin
		spawn, './' + c_prg_path + ' ' + model_name + ' ' + cfg_file
	endif else begin
		print, 'Could not found the C++ program to execute : ', c_prg_path
		print, 'Check the location of the file'
		print, 'The program will stop now'
		stop
	endelse

	; Read the outputs of c_prg_path
	output_file='model.out'
	outputs=read_outputs(output_file)
	spawn, 'rm ' + output_file ; erase the txt output file once the output is in memory

	; 'Written summary of inputs:'
	print,'Written summary of inputs:'
	print, inputs

	; plots
	e=write_on_ps_on(cfg_file)
	ymax=max(outputs[1,*])*1.1
	xr=[-2,2]*inputs.el
	plot, outputs[0,*], outputs[1,*], xtitle='Frequency', ytitle='Power',$
		yr=[0,ymax], xr=xr,/xst,/yst, background=fsc_color('white'), color=fsc_color('Black'), charsize=1.75
	for em=-inputs.el, inputs.el do begin
		plots, [inputs.fc_l - em*inputs.a1, inputs.fc_l - em*inputs.a1], [0, ymax], $
			color=fsc_color('Blue'), linestyle=2
	endfor
	if inputs.a2 ne 0 then begin
		for em=-inputs.el, inputs.el do begin  ;fc_l + m*f_s + a2_terms
			if inputs.el eq 1 then begin
				coef=(3*em*em-2)
				txt='df=a2.(3.m^2 - 2)'
			endif
			if inputs.el eq 2 then begin
				coef=(em*em-2)
				txt='df=a2.(m^2 - 2)'
			endif
			if inputs.el eq 3 then begin
				coef=(3*em*em -12)/5.
				txt='df=a2.(3.m^2 - 12)/5'
			endif
			if coef ne 0 then plots, [inputs.fc_l + (em*inputs.a1 + coef*inputs.a2), inputs.fc_l + (em*inputs.a1 + coef*inputs.a2)], [0, ymax], $
				color=fsc_color('Red'), linestyle=2
		endfor
		xyouts, 0.17, 0.9, txt, color=fsc_color('Dark Gray'), charsize=1.35, /normal
		xyouts, 0.17, 0.8, 'a2/a1='+string(100*inputs.a2/inputs.a1, format='(f5.1)')+ '%', color=fsc_color('Dark Gray'), charsize=1.35, /normal
	endif
	e=write_on_ps_off('')
	;stop
end


function read_outputs, output_file

	Nmax=2
	Kmax=200000
	out=read_Ncolumns(output_file,Nmax, Kmax, skip=skip, ref_N=0)
	data=transpose(out)
	return, data
end


function read_cfg_file, cfg_file

	cfg={range:dblarr(3), el:0, fc_l:0.0, H_l:0.0, gamma_l:0.0, a1:0.0, a2:0.0, a3:0.0, asym:0.0, inc:0.0}
	line=0
	openr, 3, cfg_file
		a=''
		while eof(3) eq 0 do begin
			readf, 3, a
			b=byte(a)
			;expected parameters and their order:  [fmin, fmax, resol], l, fc_l, H_l, gamma_l, a1, a2, a3, asym, inclination
			if a ne 35 then begin ; detect comment by searching '#' as fist caracter... ignore those lines
				line=line+1 ; the inputs must be in a specific order, we use a line number to know what parameter there
				if line eq 1 then begin ;
					readf,3,a ; read data
          			uu=strsplit(a)
          			N_uu=N_elements(uu)-1
          			r=dblarr(3)
          			for j=0,N_uu-1 do begin
          				r[j]=double(strmid(a,uu(j),uu(j+1)-uu(j)-1))
          			endfor
          				r[N_uu]=double(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))
          			cfg.range=r
				endif 
				if line eq 2 then cfg.el=long(a)
				if line eq 3 then cfg.fc_l=double(a)
				if line eq 4 then cfg.H_l=double(a)
				if line eq 5 then cfg.gamma_l=double(a)
				if line eq 6 then cfg.a1=double(a)
				if line eq 7 then cfg.a2=double(a)
				if line eq 8 then cfg.a3=double(a)
				if line eq 9 then cfg.asym=double(a)
				if line eq 10 then cfg.inc=double(a)				 
			endif
		endwhile
	close, 3

	return, cfg
end

; N: number of columns
; K: maximum number of lines
function read_Ncolumns, file,N, K, skip=skip, ref_N=ref_N

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
          		param(i,j)=double(strmid(a,uu(j),uu(j+1)-uu(j)-1))
          	endfor
			param(i, N_uu)= double(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))
		endif
		i=i+1
      endwhile

close,3
param0=param
test=where(param[*,ref_N] ne 0)

param=param[test,*]

return,param
end
