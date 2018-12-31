;---------------------------------------------------------------------------------
;Gao F, Masek J, Schwaller M, et al. On the blending of the Landsat and MODIS surface reflectance: Predicting daily Landsat surface reflectance[J]. 
;IEEE Transactions on Geoscience and Remote sensing, 2006, 44(8): 2207-2218.
;
;This code try to represent the STARFM in IDL.
;
;Input: Coarse image in t1 and t2
;Output: fine image in t1
;
;Note that:
; The size of fine image and coarse must be equal
; 
;
;parameter:
;     in the head of the main function
;
;Developed by Zhou Junxiong,email: zjxrs2018@mail.bnu.edu.cn
;
;Update history
;31/12/2018   first version
;---------------------------------------------------------------------------------

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
;-------------------------------------------------------------------------------------------------------


pro STARFM
  t0=systime(1)                  ;the initial time of program running

  ;please set the following parameters
  ;----------------------------------------------------------------------
  win=25                       ;set the half window size, if this parameter is so big, the process will be so time-consuming
  DN_min=0.0                   ;set the range of DN value of the image,If byte, 0 and 255
  DN_max=10000.0
  ;ratio=16.0                   ;set the resolution ratio, it is integer=coarse resolution/fine resolution, e.g., 480/30=16
  Fine_uncertainty = 50        ;uncertainty of fine image
  coarse_uncertainty = 50      ;uncertainty of coarse image
  resolution = 30              ;resolution of fine image
  scale_factor = 10000.0       ;define data scale factor, e.g., using MODIS reflectance products, which linearly scale refletance from 0 to 10000
  num_similar_pixel = 20       ;set number of similar pixels, too many similar pixel will blur the result
  ;------------------------------------------------------------------------

  ;open the coarse image of the first pair
  FileName1 = 'G:\FusionScale\480Data\Forest\480\T1_480'
  GetData,ImgData=coarse1,FileName = FileName1,Fid=fid1
  envi_file_query,fid1,ns=ns,nl=nl,nb=nb,dims=dims
  
  ;open the coarse image of the second pair
  FileName2 = 'G:\FusionScale\480Data\Forest\480\T2_480'
  GetData,ImgData=coarse2,FileName = FileName2,Fid=fid2

  ;open the fine img
  FileName3 = 'G:\FusionScale\480Data\Forest\T1'
  GetData,ImgData=fine,FileName = FileName3, Fid=fid3
  
  ;Assume that difference in coarse is comparable in fine
  DFusion = fine + coarse2 - coarse1
  
  ;Spectral difference between coarse and fine
  SDistance = abs(coarse1 - fine)
  SDistance_uncerntainty = sqrt((Fine_uncertainty^2+coarse_uncertainty^2))
  
  ;Temporal difference was just designed for multiple input pairs
  ;Temporal difference between coarse1 and coarse2
  ;TDistance = abs(coarse2 - coarse1)
  ;TDistance_uncerntainty = sqrt((2*coarse_uncertainty^2))
  
  ;Result of prediction
  finePre = dblarr(ns,nl,nb)
  
  ;Prediction
  for i = 0, ns-1 do begin
    for j = 0, nl-1 do begin
      
      ai=max([0,i-win])                         
      bi=min([ns-1,i+win])               
      aj=max([0,j-win])
      bj=min([nl-1,j+win])
      
      row = dblarr(bi-ai+1, bj-aj+1)
      Column = dblarr(bi-ai+1, bj-aj+1)
      for i1=0, bi-ai do begin
        Column[i1, *] = [aj:bj]
      endfor
      for j1=0, bj-aj do begin
        row[*, j1] = [ai:bi]
      endfor
      ;Spatial Distance
      SpatialDis = sqrt((row-i)^2+(column-j)^2)
      
      for k = 0, nb-1 do begin
        ;the case that do not need to use weight function
        if SDistance[i,j,k] le 10 then begin
          finePre[i,j,k] = coarse2[i,j,k]
          continue
        endif

;        if TDistance[i,j,k] le 10 then begin
;          finePre[i,j,k] = fine[i,j,k]
;          continue
;        endif
        
        ;Combined distance
        CDistance = alog(SDistance[ai:bi,aj:bj,k]+1) * SpatialDis;* alog(TDistance[ai:bi,aj:bj,k]+1) 
        
        DValue = DFusion[ai:bi, aj:bj, k]
        
        
        index = where(CDistance ne 0 and SpatialDis ne 0 and SDistance[ai:bi,aj:bj,k] lt (SDistance[i,j,k]+SDistance_uncerntainty));$
        ;and TDistance[ai:bi,aj:bj,k] lt (TDistance[i,j,k]+TDistance_uncerntainty)
        ;CDistance[index] = 1.0/CDistance[index]
        ;CDistance[index] = CDistance[index] / total(CDistance[index])
        
        ;exclude some worse neighbor pixels
        indexZero = where(CDistance eq 0 or SpatialDis eq 0 or SDistance[ai:bi,aj:bj,k] ge (SDistance[i,j,k]+SDistance_uncerntainty));$
        ;or TDistance[ai:bi,aj:bj,k] ge (TDistance[i,j,k]+TDistance_uncerntainty)
        CDistance[indexZero] = 0
        
        indcand = where(CDistance ne 0)
        CDistance = CDistance[indcand]
        order_dis=sort(CDistance)
        CDistance=CDistance[order_dis[0:num_similar_pixel-1]]           ; select the N most similar samples
        
        ;use a nomalized reverse distance as the weight function
        CDistance = 1.0/CDistance
        CDistance = CDistance / total(CDistance)
        
        DValue = DValue[indcand]
        DValue=DValue[order_dis[0:num_similar_pixel-1]]
        
        finePre[i,j,k] = total(CDistance * DValue)
      endfor
      
    endfor
  endfor

  map_info = envi_get_map_info(fid = fid3)
  OutName = 'G:\PrivateFiles\Fusion\ImageFusion\Result'
  Envi_Write_Envi_File, finePre, Out_Name = OutName, r_fid=fid_temp, ns = ns, nl = nl, nb = nb, MAP_INFO=map_info

  envi_file_mng, id = fid1, /remove
  envi_file_mng, id = fid2, /remove
  envi_file_mng, id = fid3, /remove
  ;envi_file_mng, id = fid_temp, /remove

  print, 'time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'

end