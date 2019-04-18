Function moveWinsub_TPS, coarseImg, ratio, win=win
  ;;; ratio = 480/30 = 16.0
  Mod_size = size(coarseImg)
  ns = Mod_size[1]
  nl = Mod_size[2]
  
  if KEYWORD_SET(win) eq 0 then win=5
  winsub = win-2
  nl_patch = ceil(1.0*nl/winsub)                                                               ;number of batch
  ns_patch = ceil(1.0*ns/winsub)                                                               ;
  
  coarse0_TPS = fltarr(ratio*ns, ratio*nl)
  for i=0, ns_patch-1 do begin                                                              ; iteration
    for j=0, nl_patch-1 do begin
      ai = max([0,i*winsub-1])
      aj = min([ns-1, (i+1)*winsub])
      bi = max([0,j*winsub-1])
      bj = min([nl-1, (j+1)*winsub])
      aad = 0
      bbd = 0
      if aj-ai lt 2 then begin                                                 ; 
        ai = ai-1
        aad = 1
      endif
      if bj-bi lt 2 then begin                                                 ; 
        bi = bi-1
        bbd = 1
      endif
      small_dt = coarseImg[ai:aj, bi:bj]

      ;small_dt = coarseImg[i*win:min([ns, (i+1)*win])-1, j*win:min([nl, (j+1)*win])-1]
      xgrid = [ratio/2.0, ratio]
      ygrid = [ratio/2.0, ratio]
      xout = indgen(ratio*(aj-ai+1)) + 0.5
      yout = indgen(ratio*(bj-bi+1)) + 0.5
      tmpTPS = MIN_CURVE_SURF(small_dt, xgrid=xgrid, ygrid=xgrid, /TPS, xout=xout, yout=yout)

      if (aj-ai lt win-1) or (bj-bi lt win-1) then begin                       ; 
        aii = ratio
        bii = ratio
        if aj-ai lt win-1 then begin                                           ; 
          if ai eq 0 then begin                                                ; 
            aii = 0
          endif else begin                                                     ; 
            aii = ratio
          endelse
          if aad eq 1 then begin                                               ; 
            aii = aii + ratio*1
          endif
        endif

        if bj-bi lt win-1 then begin                                           ; 
          if bi eq 0 then begin                                                ; 
            bii = 0
          endif else begin                                                     ; 
            bii = ratio
          endelse
          if bbd eq 1 then begin                                               ; 
            bii = bii + ratio*1
          endif
        endif
        tmpTPS = tmpTPS[aii:aii+ratio*min([ns, (i+1)*winsub])-i*winsub*ratio-1, bii:bii+ratio*min([nl, (j+1)*winsub])-j*winsub*ratio-1]
      endif else begin
        tmpTPS = tmpTPS[ratio:ratio*(win-1)-1, ratio:ratio*(win-1)-1]
      endelse

      coarse0_TPS[i*winsub*ratio:ratio*min([ns, (i+1)*winsub])-1, j*winsub*ratio:ratio*min([nl, (j+1)*winsub])-1] = tmpTPS
    endfor
  endfor

  return, coarse0_TPS
end