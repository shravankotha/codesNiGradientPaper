"""
    This code computes relative misorientation (misorientation of all the features wrt to one reference feature (typically the first feature))
    considering the symmetries of cubic crystals
    references: http://pajarito.materials.cmu.edu/lectures/L3-OD_symmetry-20Jan14.pdf (slide 51 for symmetry matrices)
                https://cmrl.jhu.edu/wp-content/uploads/2015/11/2006_MaterialsCharacterization_57.pdf (eq. 1 for computing misorientation)
    Input: Euler angles of features
    Output: Misorientation angles (in degrees) wrt a reference feature orientation
"""
import os
import sys
import numpy as np
import math
from alive_progress import alive_bar
from time import sleep
from computeCrystalOrientationMatrix import computeCrystalOrientationMatrix
from librarySymmetryMatricesCubic import librarySymmetryMatricesCubic
from plotHistogram import plotHistogram
import matplotlib.pyplot as plt
        
def main():

    nArguments = len(sys.argv)
    if nArguments < 3:
        text = """Two command line arguments are expected: \
                \n\t(1) input file Name (can be .ang or ASCII data from dream3D)\
                \n\t(2) output file name
               """
        raise RuntimeError(text)
        
    nameFile = str(sys.argv[1])
    out_file = str(sys.argv[2])
    #vectorRef = np.array([0,1,0])
    #arrayEulerAnglesA = np.array([0,0,0])
    #gA = computeCrystalOrientationMatrix(arrayEulerAnglesA)
     
    # find number of features in the file and also find a reference feature
    with open(nameFile, 'r') as fp:
        nLines = len(fp.readlines())
    fp.close()
    
    arrayEulerAngles = np.empty(shape=[0,2])
    with open(nameFile, 'r') as fp:
        for line in fp:
            if line[0] == '#':
                continue
            else:
                data = line.strip().split()
                if data[0].lower() == 'nan' or  data[1].lower() == 'nan' or data[2].lower() == 'nan':
                    continue
                else:
                     arrayEulerAngles_ = np.array([float(data[0]),float(data[1]),float(data[2])])
                     if (not np.any(arrayEulerAngles_)):
                        continue
                     else:
                         arrayEulerAngles = arrayEulerAngles_                         
                         break
    fp.close()
    
    file_out = open(out_file, 'w')
    if np.size(arrayEulerAngles) == 0 :
        print('***Warning: No feature with non-zero euler angles found. Taking [0 0 0] as the reference euler angles')
        arrayEulerAngles = np.array([0,0,0])
    
    gA = computeCrystalOrientationMatrix(arrayEulerAngles)                
    file_out.write('# Euler angles of reference feature : ' + str(arrayEulerAngles[0]) + '  ' + str(arrayEulerAngles[1]) + '    ' + str(arrayEulerAngles[2]) + '\n')
    listMisOri = []
    with alive_bar(nLines,title='Processing features...',bar='squares',spinner=None) as progress_bar:
        with open(nameFile, 'r') as f:
            count = 0
            for line in f:                            
                count = count + 1
                if line[0] == '#':
                    continue
                else:                    
                    data = line.strip().split()
                    if data[0].lower() == 'nan' or  data[1].lower() == 'nan' or data[2].lower() == 'nan':
                        continue
                    else:
                        arrayEulerAngles = np.array([float(data[0]),float(data[1]),float(data[2])])                        
                        gB = computeCrystalOrientationMatrix(arrayEulerAngles)
                        O = librarySymmetryMatricesCubic()
                        minMisValue = sys.float_info.max
                        for symMatrix in O:
                            tmp = (np.trace(np.dot(gA,np.matmul(np.linalg.inv(gB),symMatrix)))-1)/2
                            if abs(tmp) > 1:
                                tmp = round(tmp)
                            misOri = abs((180/math.pi)*math.acos(tmp))                                                    
                            minMisValue = min(minMisValue,misOri)
                        listMisOri.append(minMisValue)
                        file_out.write("{0:12.6f}\n".format(minMisValue))
                progress_bar()
                
    f.close()
    file_out.close()
    bins = [2*ii for ii in range(0,int(max(listMisOri)/2))]
    plotHistogram(listMisOri,bins = bins,
                  yDataType = 'fraction',
                  shouldSaveThePlot = False,
                  pathFigNameToSave = "histogram_count.png",
                  histogramBinType = 'bar',
                  binFillColor = 'black',
                  binEdgecolor = 'white', 
                  labelX = 'Angle in degrees',
                  labelY = 'Fraction of pixels',
                  scaleX = "linear",
                  scaleY = "linear",
                  fontSizeXlabel = 24,
                  fontSizeYlabel = 24,
                  fontSizeXticks = 24,
                  fontSizeYticks = 24,
                  fontName = 'Times New Roman',
                  fontWeight = "normal",
                  nPixelsPerInch = 200                  
                  )
      
if __name__ == "__main__":
    main()    