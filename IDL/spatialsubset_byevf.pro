;+
; 
; DESCRIPTION:
; This pro can subset raster by evf
; Can't deal with emtypyhole
; 
;SYNTAX
;   spatialsubset_byevf, data_fid, evfname, out_name
;
;INPUT PARAMETERS:
;   data_fid: FID of the raster image
;   evfname: the file path of evf file
;   out_name: Filename of result
;
;AUTHOR: ZJX @ Beijing Nomral University (2018/12/08) [zjxrs2014@163.com]
; 
;-


pro spatialsubset_byevf,data_fid,evfname,out_name
  compile_opt idl2
  envi_file_query,data_fid,bnames=bnames,ns=ns,nl=nl,nb=nb
  
  ;open evf
  evf_id = envi_evf_open(evfname)
  ;get info of the evf
  envi_evf_info,evf_id,num_recs=num_recs,data_type=data_type,projection=projection,layer_name=layer_name
  roi_ids = lonarr(num_recs)
  ;
  ;oproj = envi_get_projection(fid=data_fid)
  
  for i=0,num_recs-1 do begin
    record = envi_evf_read_record(evf_id,i)
    ;According to image geo and evf geo
    ;envi_convert_projection_coordinates,record[0,*],record[1,*],projection,oxmap,oymap,oproj
    envi_convert_file_coordinates,data_fid,xmap,ymap,record[0,*],record[1,*]
    
    ;create ROI
    roi_id = envi_create_roi(color=4,ns=ns,nl=nl)
    envi_define_roi,roi_id,/polygon,xpts=reform(xmap),ypts=reform(ymap)
    roi_ids[i] = roi_id
    
    ;Spatial extent
    if i eq 0 then begin
      xmin = round(min(xmap,max=xmax))
      ymin = round(min(ymap,max=ymax))
    endif else begin
       xmin = xmin < round(min(xmap))
       xmax = xmax > round(max(xmap))
       ymin = ymin < round(min(ymap))
       ymax = ymax > round(max(ymap)) 
    endelse
  endfor
  xmin = xmin > 0
  xmax = xmax < ns-1
  ymin = ymin > 0
  ymax = ymax < nl-1
  
  ;Create Mask
  envi_mask_doit,and_or=1,/in_memory,roi_ids=roi_ids,ns=ns,nl=nl,/inside,r_fid=m_fid
  out_dims = [-1,xmin,xmax,ymin,ymax]
  
  ;Subset
  envi_mask_apply_doit,fid=data_fid,pos=indgen(nb),dims=out_dims,m_fid=m_fid,m_pos=[0],value=0,out_bname=bnames,$
    in_memory=0,out_name=out_name,r_fid=r_fid
  ;掩膜文件ID移除
  envi_file_mng,id=m_fid,/remove
end