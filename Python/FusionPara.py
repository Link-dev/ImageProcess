
from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt


def calcAB(x, y):
    
    size = len(x)
    a = 0.0
    b = 0.0

    a = (size * sum((x * y)) - sum(x) * sum(y)) / (size * sum(x ** 2) - (sum(x)) ** 2)
    b = np.mean(y) - a * np.mean(x)
    return a, b


#FSDAF
#'D://FusionScale//480Data//Coleambally//', 'D://FusionScale//480Data//Forest//', \
#           'D://FusionScale//480Data//Gwydir//', 'D://FusionScale//480Data//Heterogeneous//', \

#FileName = ['D://FusionScale//480Data//Coleambally//', 'D://FusionScale//480Data//Forest//', \
#            'D://FusionScale//480Data//Gwydir//', 'D://FusionScale//480Data//Heterogeneous//']

# 
FileName = ['D://FusionScale//480Data//Coleambally//', 'D://FusionScale//480Data//Forest//',\
            'D://FusionScale//480Data//Gwydir//', 'D://FusionScale//480Data//Heterogeneous//', \
            'D://FusionScale//480Data//NewSite//Site1//', 'D://FusionScale//480Data//NewSite//Site2//', \
            'D://FusionScale//480Data//NewSite//Site3//', 'D://FusionScale//480Data//NewSite//Site4//']



NumFile = len(FileName)

scale = [2, 4, 8, 16, 32, 40]
NumScale = len(scale)

#波段数
nb = 1
ImgXSize = 480
ImgYSize = 480

RMSE = np.zeros(NumFile * NumScale * nb)
F_I = np.zeros(NumFile * NumScale * nb)

Index = 0

preN = 2

for iFile in range(0, NumFile):
    File1 = FileName[iFile] + 'T1'
    #File2 = FileName[iFile] + 'T2'
    if (iFile < 4):
        preN = 2
    else:
        preN = 3
    
    File2 = FileName[iFile] + 'T' + str(preN)
    
    Dataset1 = gdal.Open(File1)
    Dataset2 = gdal.Open(File2)

    for jScale in range(0, NumScale):
        File3 = FileName[iFile]  + str(scale[jScale] * 30) + '//T1_' + str(scale[jScale] * 30)
        File4 = FileName[iFile] + str(scale[jScale] * 30) + '//T' + str(preN) + '_' + str(scale[jScale] * 30)
        #File5 = FileName[iFile] + 'STARFM//Result//pre_' + str(scale[jScale] * 30)
        File5 = FileName[iFile] + str(scale[jScale] * 30) + '//T' + str(preN) + '_' + str(scale[jScale] * 30) + '_FSDAF'
        Dataset3 = gdal.Open(File3)
        Dataset4 = gdal.Open(File4)
        Dataset5 = gdal.Open(File5)

        RMSE_S = np.zeros(nb)
        Fstd_S = np.zeros(nb)
        F_C = np.zeros(nb)
        for iBand in range(0, nb):
            band = Dataset1.GetRasterBand(iBand + 1)
            data1 = band.ReadAsArray(0, 0, ImgXSize, ImgYSize)/10000.0

            band = Dataset2.GetRasterBand(iBand + 1)
            data2 = band.ReadAsArray(0, 0, ImgXSize, ImgYSize) / 10000.0

            band = Dataset3.GetRasterBand(iBand + 1)
            data3 = band.ReadAsArray(0, 0, ImgXSize, ImgYSize) / 10000.0

            band = Dataset4.GetRasterBand(iBand + 1)
            data4 = band.ReadAsArray(0, 0, ImgXSize, ImgYSize) / 10000.0

            band = Dataset5.GetRasterBand(iBand + 1)
            data5 = band.ReadAsArray(0, 0, ImgXSize, ImgYSize) / 10000.0

            data = (data5 - data2) ** 2
            RMSE_S[iBand] = np.sqrt(np.mean(data))

            imgSize = int(ImgXSize/scale[jScale]*ImgYSize/scale[jScale])
            CoarseXSize = int(ImgXSize/scale[jScale])
            CoarseYSize = int(ImgYSize/scale[jScale])
            #RMSE_C = np.zeros(imgSize)
            #Fstd_C = np.zeros(imgSize)
            Fvar_C = np.zeros(imgSize)
            
            for j in range(0, CoarseXSize, 1):
                for k in range(0, CoarseYSize, 1):
                    tmp1 = data1[[x for x in range(j, j + scale[jScale])]]
                    tmp1 = tmp1[:, [x for x in range(k, k + scale[jScale])]]

                    tmp2 = data2[[x for x in range(j, j + scale[jScale])]]
                    tmp2 = tmp2[:, [x for x in range(k, k + scale[jScale])]]

                    tmp5 = data5[[x for x in range(j, j + scale[jScale])]]
                    tmp5 = tmp5[:, [x for x in range(k, k + scale[jScale])]]

                    tmp = (tmp1 - tmp2)
                    Fvar_C[j*CoarseYSize + k] = np.var(tmp)
                    #Fstd_C[j*CoarseYSize + k] = np.std(tmp)
                    
                    
                    #RMSE_C[j*CoarseYSize + k] = np.sqrt(sum(sum((tmp5 - tmp2) ** 2))/imgSize)

            #RMSE与变化量标准差的关系，
            """
            a, b = calcAB(Fstd_C, RMSE_C)
            print("y = %10.5fx + %10.5f" % (a, b))
            maxL = np.max(Fstd_C)
            minL = np.min(Fstd_C)
            x1 = np.linspace(minL, maxL)
            y1 = a * x1 + b
            plt.plot(x1, y1, 'r')
            plt.scatter(Fstd_C, RMSE_C)
            R2 = 1 - sum((RMSE_C - (a * Fstd_C + b)) ** 2) / sum((RMSE_C - np.mean(RMSE_C)) ** 2)
            print(R2)
            plt.show()
            #RMSE = np.sqrt(sum((RMSE_C - (a * Fstd_C + b)) ** 2) / len(Fstd_C))
            """
                    
            Fstd_S[iBand] = np.sqrt(np.mean(Fvar_C))
            F_C[iBand] = np.mean(abs(data4 - data3))
            
            Index = Index + 1

        RMSE[Index-nb:Index] = RMSE_S
        F_I[Index-nb:Index] = np.log(F_C * Fstd_S)

a, b = calcAB(F_I, RMSE)
print("y = %10.5fx + %10.5f" % (a, b))
maxL = np.max(F_I)
minL = np.min(F_I)
x1 = np.linspace(minL, maxL)
y1 = a * x1 + b
plt.plot(x1, y1, 'r')
plt.scatter(F_I, RMSE)
R2 = 1 - sum((RMSE - (a * F_I + b)) ** 2) / sum((RMSE - np.mean(RMSE)) ** 2)
print(R2)
plt.show()
