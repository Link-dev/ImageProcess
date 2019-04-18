
; -------------------------------------------------------------------------------------------------------------------
; universal function: function for open the file
; -------------------------------------------------------------------------------------------------------------------
Pro GetData,ImgData = ImgData,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,$
  FileName = FileName,Map_info = map_Info, Fid = Fid, dims = dims
  Envi_Open_File,FileName,R_Fid = Fid
  Envi_File_Query,Fid,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type
  map_info = envi_get_map_info(fid=Fid)
  dims = [-1,0,ns - 1 ,0,nl - 1]
  case Data_Type Of
    1:ImgData = BytArr(ns,nl,nb)    ;  BYTE  Byte
    2:ImgData = IntArr(ns,nl,nb)    ;  INT  Integer
    3:ImgData = LonArr(ns,nl,nb)    ;  LONG  Longword integer
    4:ImgData = FltArr(ns,nl,nb)    ;  FLOAT  Floating point
    5:ImgData = DblArr(ns,nl,nb)    ;  DOUBLE  Double-precision floating
    6:ImgData = COMPLEXARR(ns,nl,nb); complex, single-precision, floating-point
    9:ImgData = DCOMPLEXARR(ns,nl,nb);complex, double-precision, floating-point
    12:ImgData = UINTARR(ns,nl,nb)   ; unsigned integer vector or array
    13:ImgData = ULONARR(ns,nl,nb)   ;  unsigned longword integer vector or array
    14:ImgData = LON64ARR(ns,nl,nb)   ;a 64-bit integer vector or array
    15:ImgData = ULON64ARR(ns,nl,nb)   ;an unsigned 64-bit integer vector or array
  EndCase
  For i = 0,nb-1 Do Begin
    Dt = Envi_Get_Data(Fid = Fid,dims = dims,pos=i)
    ImgData[*,*,i] = Dt[*,*]
  EndFor
End
;------------------------------------------------------------------------------------------------------

function tps_parallel, data, scale, win, threads, path
  w = win/2
  ;---------parameter setting----------------------------------
  threads=threads ; !!number of thread
  Function_position=path + 'moveWinsub_TPS.pro' ;!! Function of child thread
  ;------------------------------------------------------------

  ImgData_size=size(data)
  ns=ImgData_size[1]
  nl=ImgData_size[2]

  result=fltarr(ns*scale,nl*scale)

  p=objarr(threads) ; objects array for threads
  task_num=ceil(nl/(win-4)) - 1; tasks number
  new_task=0

  ;------------------thread initialization-------------------------------
  FOR ti=0,threads-1 do begin

    ; Establish error handler. When errors occur, the index of the

    ; error is returned in the variable Error_status:

    CATCH, Error_status

    ;This statement begins the error handler:

    IF Error_status NE 0 THEN BEGIN

      PRINT, 'Error index: ', Error_status

      PRINT, 'Error message: ', !ERROR_STATE.MSG

      ; Handle the error by extending A:

      wait, 60

      p[ti]=OBJ_NEW('IDL_IDLBridge')

      CATCH, /CANCEL

    ENDIF

    p[ti]=OBJ_NEW('IDL_IDLBridge')
    p[ti]->Execute,".compile "+"'"+Function_position+"'"
    p[ti]->Execute,".compile "+"'"+Function_position+"'"

    Data_p=Data[*,new_task*(win-4):new_task*(win-4)+win-1]
    p[ti]->Setvar,"Data_p",Data_p
    p[ti]->Setvar,"scale",scale
    p[ti]->Setvar,"j",new_task
    new_task++
    p[ti]->Execute,"Fine=moveWinsub_TPS(Data_p, scale)",/nowait
  ENDFOR

  signal=intarr(threads); threads' status managernew_task=0
  ;print,'Successfully initialized and waiting for signal'


  ;---------------------Running---------------------------------------------
  While (1 gt 0) do begin
    ;check the slave，1 is busy，0 is idle，2 is end
    for tj=0,threads-1 do signal[tj]=p[tj]->Status()
    ;if all slaves are idle, then break
    if (new_task eq task_num)and(total(signal)eq 0) then break;all tasks are done
    pos=where(signal eq 0 or signal eq 2, pos_count);idle threads
    
    ;assign new task to the idle slaves
    FOR ti=0,pos_count-1 do begin
      thread_idle=pos[ti]
      ;
      tj=thread_idle
      j=p[tj]->Getvar('j')

      result1 = p[tj]->Getvar('Fine')

      ;deal with the edge pixels
      if (j ne 0 && j ne task_num-1) then begin
        result[*, (j*(win-4)+2)*scale:((j*(win-4)+2*w-1)*scale-1)]=result1[*, 2*scale: (win-2)*scale-1]
      endif else if j eq 0 then begin
        result[*, 0:((j*(win-4)+2*w-1)*scale-1)] = result1[*, 0:(win-2)*scale-1]
      endif else if j eq task_num-1 then begin
        Data_size=size(result1)
        Data_nl=Data_size[2]
        result[*, (j*(win-4)+2)*scale:nl*scale-1] = result1[*, 2*scale:Data_nl-1]
      endif

      if (new_task eq task_num) then continue;no tasks left

      if new_task ne task_num-1 then begin
        Data_p=Data[*,new_task*(win-4):new_task*(win-4)+win-1]
      endif else begin
        Data_p=Data[*,new_task*(win-4):nl-1]
      endelse

      p[tj]->Setvar,"Data_p",Data_p
      p[tj]->Setvar,"scale",scale
      p[tj]->Setvar,"j",new_task
      new_task++
      p[tj]->Execute,"Fine=moveWinsub_TPS(Data_p, scale)",/nowait

      ;print,'Moves to task:'+string(new_task)
      ;print,'thread:'+string(tj)

      ;print,'Assigned/Unassigned task: '+string(new_task,format='(i)')+'/'+string(task_num,format='(i)')
      ;print,'*******************************************************************'
      ;print,''
    ENDFOR
  ENDwhile

  ;release all the slaves
  FOR ti=0,threads-1 do begin
    p[ti]->cleanup
  ENDFOR

  ;print,'Computation is completed!'
  return,result
end

;-------------------------------------------------------------------
;
;                       main program
;-------------------------------------------------------------------
;
;
;This code is for tps interpolation by parallel, but it may not be good in small image
;And this code perform bad due to block effects
;Just a sample for IDL parallel
;date: 20190418

pro  TPS_project_Final
  envi, /restore_base_save_files
  envi_batch_init

  t0=systime(1)             ; the initial time of program running

  ; please set the following parameters
  ;----------------------------------------------------------------------
  ratio=8                   ; a MODIS pixel(240m) consists of 8*8 Landsat pixels
  win = 7                   ; window size for interpolation
  threads = 8               ; number of thread
  ;------------------------------------------------------------------------
  
  path = 'G:\final\'
  
  ;open the fine image
  FileName1 = path + 'L7_2001_10_07_flaash_CIA_patch_img_upscale'
  GetData,ImgData=Coarse,FileName = FileName1,Fid=fid
  envi_file_query,fid,ns=ns_coarse,nl=nl_coarse,nb=nb_coarse,dims=dims, wl = wavelength
  
  map_info = envi_get_map_info(fid=fid)
  stackBandName = strarr(nb_coarse)
  
  
  TPS = fltarr(ns_coarse*ratio, nl_coarse*ratio, nb_coarse)
  
  for i = 0, nb_coarse-1 do begin
    band = Coarse[*,*,i]
    TPS[*,*, i] = tps_parallel(band, ratio, win, threads, path)

    stackBandName[i] = 'TPS' + strtrim(i+1,1)
  endfor
  Envi_Write_Envi_File, TPS, Out_Name = path + '\TPS_Prediction', r_fid=fid1, map_info = map_info, wl = wavelength
  
  envi_file_mng, id=fid, /remove
  ;envi_file_mng, id=fid1, /remove
  
  print, floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'
  
end
