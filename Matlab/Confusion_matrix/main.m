clc;clear;
%%数据路径
path1 = 'G:\PrivateFiles\shuju\分类结果\';
path2 = 'G:\PrivateFiles\shuju\真值\';

%数据名称
classMapName = '16GLCM.tif'; %分类结果
TrueMapName = '16g.tif'; %真值

classMapName = [path1 classMapName];
TrueMapName = [path2 TrueMapName];

classMap = imread(classMapName);
TrueMap = imread(TrueMapName);

classMap = double(classMap(:));
TrueMap = double(TrueMap(:));

confusion_matrix(TrueMap, classMap);