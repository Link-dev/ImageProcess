
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

pro bandStackingFun, scale_factor = scale_factor, folder = folder


ratio = scale_factor

FileName1 = folder + 'T1'
GetData,ImgData=fine0,FileName = FileName1, Fid=fid0
 envi_file_query,fid0,ns=ns_fine,nl=nl_fine,nb=nb,dims=dims_fine

fine = fltarr(ns_fine, nl_fine, nb)
stackBandName = strarr(nb)

folder = folder + 'STARFM\Result\'
for iband = 0, nb-1 do begin
  
  FileName1 = folder + 'pre_' + strtrim(ratio * 30,1) + 'band' + strtrim(iband+1,1) + '.bin'
  GetData,ImgData=fine1,FileName = FileName1, Fid=fid1

  stackBandName[iband] = strtrim(iband+1,1)
  
  tmp1 = fine1
  
  inot = where(fine1[*,*] lt 0, num_back)
  
  if (num_back ne 0) then begin
    tmp2 = fine0[*,*,iband]
    tmp1[inot] = tmp2[inot]
  endif
  
  fine[*,*,iband] = tmp1
  
  envi_file_mng, id = fid1, /remove

endfor



;输出融合结果
stackName = folder + 'pre_' + strtrim(ratio * 30,1)

Envi_Write_Envi_File, fine, Out_Name = stackName, r_fid=fid_temp, ns = ns_fine, nl = nl_fine, nb = nb, out_dt = Data_Type

end

pro bandstacking
  scale = [48, 60, 80, 96, 120]
  NumScale = n_elements(scale)

  folder=['D:\FusionScale\480Data\Coleambally\', 'D:\FusionScale\480Data\Forest\', 'D:\FusionScale\480Data\Gwydir\', 'D:\FusionScale\480Data\Heterogeneous\']
  ;folder = ['D:\FusionScale\480Data\NewSite\Site1\', 'D:\FusionScale\480Data\NewSite\Site2\', 'D:\FusionScale\480Data\NewSite\Site3\', 'D:\FusionScale\480Data\NewSite\Site4\']
  for i = 0, 3, 1 do begin
    for j = 0, NumScale-1, 1 do begin
      bandStackingFun,scale_factor = scale[j], folder = folder[i]
    endfor
  endfor
end