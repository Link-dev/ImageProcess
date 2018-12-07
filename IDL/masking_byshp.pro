;+
;                        masking_byshp.pro
;DESCRIPTION:
; This pro is the example of subset by shp/evf
;
;INPUT PARAMETERS:
;   filename: the file path of raster image
;   shp_file: the file path of shp file
;   outName: Filename of result
;   
;AUTHOR: ZJX @ Beijing Nomral University (2018/12) [zjxrs2014@163.com]
;        
;-

pro masking_byshp
  compile_opt idl2
  ;envi,/restore_base_save_files
  envi_batch_init
  
  ;open the origin image
  filename = 'G:\code\FileRecv\relative_eff_light_2001_2012.tif'
  Envi_Open_File,filename,r_fid=fid
  
  ;shapefile
  shp_file = 'G:\code\FileRecv\2000_shp\jjj_maincities_2000.shp'
  
  ;the fileFolder of output
  outName = 'G:\code\FileRecv\SubsetImg.tif'
  
  ;ROI_file_name = 'Q:\code\FileRecv\roi.roi'
  ;ROIhole_file_name = 'Q:\code\FileRecv\roiHole.roi'
  
  ;convert shp to roi
  spatialsubset_byshp, fid, shp_file, outName
  

end