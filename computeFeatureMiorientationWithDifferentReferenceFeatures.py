"""
    This code does the following: (Essentially computing the misorientation of one random pixel wrt to large number of other sampled pixels and finally computing the statistics of this data)
    (1) computes relative misorientation (misorientation of nSamplesInEachTrial number of features wrt to a RANDOMLY CHOSEN reference feature)
    (2) The misorientation angle corresponding to requiredMisorientationPercentile is computed (This is basically the angle below which the misorientations of requiredMisorientationPercentile percentage of features lie)
    (3) This is repeated for a number of trials and histogram of misOriPercentile is plotted    
    Input: Euler angles of features
    Output: misorientation angle corresonding to the given percentile in each trial is written to a file
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
import random
        
def main():

    nArguments = len(sys.argv)
    if nArguments < 6:
        text = """Five command line arguments are expected: \
                \n\t(1) input file Name (Euler angles file)\
                \n\t(2) number of header lines \
                \n\t(3) required Misorientation Percentile (0-100) \
                \n\t(4) number of trials to compute the percentiles \
                \n\t(5) number of samples in each trial \
                \n\t(6) output file to write percentiles
               """
        raise RuntimeError(text)
        
    nameFile = str(sys.argv[1])
    nHeaderLines = int(sys.argv[2])
    requiredMisorientationPercentile = float(sys.argv[3])
    nTrials = int(sys.argv[4])
    nSamplesInEachTrial = int(sys.argv[5])
    out_file = str(sys.argv[6])
    tolerance = 1E-6
    maxNumberOfAttempsToFindRefPixel = 20

    dataFile = np.genfromtxt(nameFile, skip_header = nHeaderLines, dtype = 'str')
    nRowsEulerAngles = np.shape(dataFile)[0] - 1
    file_out = open(out_file, 'w') 
    file_out.write('trialNumber     misOriAngleForTheGivenPercentile \n')    

    listMisOri, listPercentiles = [[] for ii in range(0,2)]
    with alive_bar(nTrials,title='Processing trials...',bar='squares',spinner=None) as progress_bar:    
        for iTrial in range(0,nTrials):
            # randomly chose a reference pixel
            idRowRef = random.randint(0,nRowsEulerAngles)
            isNotFoundAValidRefPixel = True
            nAttemps = 0
            while isNotFoundAValidRefPixel and nAttemps < maxNumberOfAttempsToFindRefPixel:
                nAttemps = nAttemps + 1
                if dataFile[idRowRef,0].lower() == 'nan' or  dataFile[idRowRef,1].lower() == 'nan' or dataFile[idRowRef,2].lower() == 'nan' or  \
                (float(dataFile[idRowRef,0]) < tolerance and  float(dataFile[idRowRef,1]) < tolerance and float(dataFile[idRowRef,2]) < tolerance):
                    isNotFoundAValidRefPixel = True
                else:
                    arrayEulerAngles = np.array([float(dataFile[idRowRef,0]),float(dataFile[idRowRef,1]),float(dataFile[idRowRef,2])])
                    isNotFoundAValidRefPixel = False
            if isNotFoundAValidRefPixel:
                raise RuntimeError("No valid reference pixel is found in the maximum number of attempts")
                
            gA = computeCrystalOrientationMatrix(arrayEulerAngles)    
    
            listMisOri = []
            for jj in range(0,nSamplesInEachTrial):
                ii = random.randint(0,nRowsEulerAngles)
                if dataFile[ii,0].lower() == 'nan' or  dataFile[ii,1].lower() == 'nan' or dataFile[ii,2].lower() == 'nan':
                    continue
                else:
                    arrayEulerAngles = np.array([float(dataFile[ii,0]),float(dataFile[ii,1]),float(dataFile[ii,2])])
                    gB = computeCrystalOrientationMatrix(arrayEulerAngles)
                    O = librarySymmetryMatricesCubic()
                    minMisValue = sys.float_info.max
                    for symMatrix in O:
                        misOri = abs((180/math.pi)*math.acos((np.trace(np.dot(gA,np.matmul(np.linalg.inv(gB),symMatrix)))-1)/2))
                        minMisValue = min(minMisValue,misOri)
                    listMisOri.append(minMisValue)
            misOriPercentile = np.percentile(listMisOri, requiredMisorientationPercentile)        
            listPercentiles.append(misOriPercentile)
            file_out.write("{0:10d}{1:12.6f}\n".format(iTrial,misOriPercentile))
            file_out.flush()
            progress_bar()
                
    file_out.close()
    bins = [2*ii for ii in range(0,int(max(listPercentiles)/2))]
    plotHistogram(listPercentiles,bins = bins,
                  yDataType = 'fraction',
                  shouldSaveThePlot = False,
                  pathFigNameToSave = "histogram_count.png",
                  histogramBinType = 'bar',
                  binFillColor = 'black',
                  binEdgecolor = 'white', 
                  labelX = 'Angle in degrees',
                  labelY = 'Fraction of trials',
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