clc;clear;
%%����·��
path1 = 'G:\PrivateFiles\shuju\������\';
path2 = 'G:\PrivateFiles\shuju\��ֵ\';

%��������
classMapName = '16GLCM.tif'; %������
TrueMapName = '16g.tif'; %��ֵ

classMapName = [path1 classMapName];
TrueMapName = [path2 TrueMapName];

classMap = imread(classMapName);
TrueMap = imread(TrueMapName);

classMap = double(classMap(:));
TrueMap = double(TrueMap(:));

confusion_matrix(TrueMap, classMap);