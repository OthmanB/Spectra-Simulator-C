; N: number of columns
; K: maximum number of lines
function read_txt, file,Ncols, Nlines, skip

;skip=0

if n_elements(skip) eq 0 then skip=1 ; defaut we skip one line only

openr, 3, file

	param=dblarr(Nlines,Ncols)
	a=''
	i=0d
	ii=0d
      while EOF(3) ne 1 do begin
      	if i lt skip then readf,3,format='(q)'
        if i ge skip then begin
          	readf,3,a ; read data
          	uu=strsplit(a)
          	N_uu=N_elements(uu)-1
          	for j=0,N_uu-1 do begin
          		param(ii,j)=float(strmid(a,uu(j),uu(j+1)-uu(j)-1))
          	endfor
			param(ii, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))
			ii=ii+1
		endif
		i=i+1
      endwhile

close,3
param=transpose(param[0:ii-1,*])

return,param
end

