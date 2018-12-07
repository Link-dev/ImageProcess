;+
;
; DESCRIPTION:
; This pro can subset raster by shapefile
; Can deal with emtypyhole, but need to be revised
;
;SYNTAX
;   spatialsubset_byevf, data_fid, evfname, out_name
;
;INPUT PARAMETERS:
;   fid: FID of the raster image
;   shp_file: the file path of shp file
;   out_name: Filename of result
;
;AUTHOR: ZJX @ Beijing Nomral University (2018/12/08) [zjxrs2014@163.com]
;
;-

;This function can deal with blank hole, but I don't know why it can't work now
Function clockwise,polygon_verts
  ;---------------------------------------------------------------
  ;DESCRIPTION:
  ;   This function is to determine the polygon verts is whether clockwise or not.
  ;
  ;SYNTAX
  ;   Result = clockwise (polygon_verts)
  ;
  ;RETURN VALUE: 1 (the polygon verts is clockwise )or 0 (the polygon verts is not clockwise)
  ;
  ;INPUT PARAMETER:
  ;   polygon_verts: 2-dimension matrix (2*verts) denoting the (X,Y) recodes of polygon verts
  ;---------------------------------------------------------------

  vert_num=n_elements(polygon_verts)/2
  Y_value=polygon_verts[1,*]
  i_vert=where (Y_value eq max(Y_value), v_num)
  a=0
  for iv=0,v_num-1 do begin
    if (i_vert[iv] ne 0) and (i_vert[iv] ne vert_num-1) then begin
      kk=(polygon_verts[0,i_vert[iv]]-polygon_verts[0,i_vert[iv]-1])*(polygon_verts[1,i_vert[iv]+1]-polygon_verts[1,i_vert[iv]])- $
        (polygon_verts[1,i_vert[iv]]-polygon_verts[1,i_vert[iv]-1])*(polygon_verts[0,i_vert[iv]+1]-polygon_verts[0,i_vert[iv]])
    endif
    if (i_vert[iv] eq 0) or (i_vert[iv] eq vert_num-1) then begin
      kk=(polygon_verts[0,0]-polygon_verts[0,vert_num-2])*(polygon_verts[1,1]-polygon_verts[1,0])- $
        (polygon_verts[1,0]-polygon_verts[1,vert_num-2])*(polygon_verts[0,1]-polygon_verts[0,0])
    endif
    a=a+(kk le 0)
  endfor

  clock= a gt 0
  return, clock

END

pro spatialsubset_byshp, fid, shp_file, outName

  IF (keyword_set(tagname_num) eq 0) THEN tagname_num=0L
  envi_file_query, fid, fname=fname, ns=ns, nl=nl, nb=nb, bnames=bnames,data_type=data_type
  map_info = envi_get_map_info(fid=fid)
  oproj = ENVI_GET_PROJECTION(FID = fid)

  ;Read the projection information
  shp_file_len=strlen(shp_file)
  shp_basename=file_basename(shp_file,'.shp')
  shp_name=STRMID(shp_file,0,shp_file_len-4)
  projstr=shp_name+'.prj'
  openr,lun,projstr,/GET_LUN
  shpPrjstr=''
  readf,lun,shpPrjstr
  free_lun,lun
  shapeproject=envi_proj_create(type=42,pe_coord_sys_str=shpPrjstr)
  if shapeproject.type eq 0 then shapeproject=envi_proj_create(type=1,pe_coord_sys_str=shpPrjstr)

  oshp=Obj_New('IDLffShape',shp_file)
  oshp->getproperty,n_entities=n_ent,attribute_info=attr_info,$
    n_attributes=n_attr,Entity_type=ent_type
  attr_tag=attr_info.name
  print,attr_tag;attr_struct={}
  attr=oshp->GetAttributes(/All)

  ROI_IDS=intarr(n_ent)
  ROIhole_IDS=intarr(n_ent)
  for i=0,n_ent-1 do begin
    ent=oshp->GetEntity(i)
    N_parts=ent.N_PARTS
    vert=*(ent.vertices)

    name=STRTRIM(string(attr[i].(tagname_num)),2); NEED TO BE CHANGED! the name of ROI

    ROI_ID=ENVI_CREATE_ROI(color=i*2+2,name=name, ns = ns ,  nl = nl)
    ROI_ID_hole=ENVI_CREATE_ROI(color=i*2+3,name=name+'_hole', ns = ns ,  nl = nl)
    if N_parts le 1 then begin
      
      ;if projection exists, use this
      envi_convert_projection_coordinates, vert[0,*], vert[1,*], shapeproject, oxmap, oymap, oproj
      ENVI_CONVERT_FILE_COORDINATES,fid,xpts,ypts, oxmap, oymap
      
      ;ENVI_CONVERT_FILE_COORDINATES,fid,xmap, ymap, vert[0,*], vert[1,*]
      
      ENVI_DEFINE_ROI, roi_id, /polygon, xpts=reform(xpts), ypts=reform(ypts)
      
    endif else begin
      N_VERTICES=ent.N_VERTICES
      parts=*(ent.parts)

      if N_parts ge 2 then begin
        ;parts_ptr=[0,parts[1]]
        for part_i=0, N_parts-1 do begin
          if  part_i ne N_parts-1 then begin
            polygon_verts=vert[*,parts[part_i]:(parts[part_i+1]-1)]
            envi_convert_projection_coordinates, polygon_verts[0,*], polygon_verts[1,*], shapeproject, oxmap, oymap, oproj
            ENVI_CONVERT_FILE_COORDINATES,fid,xpts,ypts,oxmap,oymap

            if clockwise(polygon_verts) eq 1 then ENVI_DEFINE_ROI, roi_id, /polygon, xpts=reform(xpts), ypts=reform(ypts)
            if clockwise(polygon_verts) eq 0 then ENVI_DEFINE_ROI, roi_id_hole, /polygon, xpts=reform(xpts), ypts=reform(ypts)

          endif

          if part_i eq N_parts-1 then begin
            polygon_verts=vert[*,parts[part_i]:(N_VERTICES-1)]
            envi_convert_projection_coordinates, polygon_verts[0,*], polygon_verts[1,*], shapeproject, oxmap, oymap, oproj
            ENVI_CONVERT_FILE_COORDINATES,fid,xpts,ypts, oxmap, oymap
            if clockwise(polygon_verts) eq 1 then ENVI_DEFINE_ROI, roi_id, /polygon, xpts=reform(xpts), ypts=reform(ypts)
            if clockwise(polygon_verts) eq 0 then ENVI_DEFINE_ROI, roi_id_hole, /polygon, xpts=reform(xpts), ypts=reform(ypts)

          endif
        endfor
      endif
      
    endelse
    roi_ids[i]=roi_id
    ROIhole_IDS[i]=roi_id_hole
    
    ;for clipping
    if i eq 0 then begin
      xmin = round(min(xpts,max=xmax))
      ymin = round(min(ypts,max=ymax))
    endif else begin
      xmin = xmin < round(min(xpts))
      xmax = xmax > round(max(xpts))
      ymin = ymin < round(min(ypts))
      ymax = ymax > round(max(ypts))
    endelse
    
  endfor
  xmin = xmin > 0
  xmax = xmax < ns-1
  ymin = ymin > 0
  ymax = ymax < nl-1
  ;创建掩膜，裁剪后掩
  envi_mask_doit,and_or=1,/in_memory,roi_ids=roi_ids,ns=ns,nl=nl,/inside,r_fid=m_fid
  ;envi_mask_doit,and_or=1,/in_memory,roi_ids=ROIhole_IDS,ns=ns,nl=nl,/inside,r_fid=m_fid1
  ;out_dims = [-1, 0, ns-1, 0, nl-1]
  
  out_dims1 = [-1,xmin,xmax,ymin,ymax]
  envi_mask_apply_doit,fid=fid, pos = [0], dims=out_dims, m_fid=m_fid,m_pos=[0],value=0,$
    in_memory=0,out_name=outName,r_fid=r_fid
  
  
  ;envi_mask_apply_doit,fid=r_fid, pos = [0], dims=out_dims1,m_fid=m_fid1,m_pos=[0],value=1,$
  ;  in_memory=0,out_name=outName,r_fid=r_fid1
  ;掩膜文件ID移除
  ;envi_file_mng,id=m_fid,/remove
  
  ;ENVI_SAVE_ROIS, ROI_file_name, ROI_IDS
  ;ENVI_SAVE_ROIS, ROIhole_file_name, ROIhole_IDS
end