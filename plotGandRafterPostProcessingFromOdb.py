''' This codes plots various quantities as given below for different layers in the deposit:
    1. Thermal gradient
    2. Solidification velocity
    3. Cooling rate
    4. PDAS etc.
'''

import sys
import os
import statistics as stats
import random as rand
import time
from plotScatter4D import plotScatter4D
from plotScatter2D import plotScatter2D, plotScatter2DwithErrorBars
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import math
from functionComputePDAS import functionComputePDAS_NiSuperAlloys as computePDAS_NiSuperAlloys


def main():    
    nArguments = len(sys.argv)
    if nArguments != 3:
        text = """Two command line arguments are expected: \
                \n\t(1) file name \
                \n\t(2) gradingDistribution (625centerTo738/738centerTo625) """
        raise RuntimeError(text)

    nameFile = str(sys.argv[1])
    gradingDistribution = str(sys.argv[2])
    pwd = os.getcwd()
    
    imgFormat                   = ".tiff" # raster -- .png,.pdf,.jpg,.tiff; vector -- .svg
    nPixelsPerInch              = 1000
    fontName                    = "Times New Roman"
    legendFont                  = fm.FontProperties(family=fontName,size=24, weight='normal')
    fontStyleForLabels          = {'fontname':fontName, 'fontsize':24,'weight':'normal'}
    fontDictPlotText            = {'fontname':fontName,'fontsize':24,'weight':'normal'}
    fontSizeXticks              = 24
    fontSizeYticks              = 24    
    fontWeightXticks            = 'normal'
    fontWeightYticks            = 'normal'
    lineWidth                   = 2
    plotScatterPlotsIn3D        = False
    plotRegionsSeparately       = False
    plotLayersSeparately        = False
    plotRadialDataSeparately    = False
    
    nPointsToPlot               = 10000
    multiplierThermalGradient   = 1000
    idColumnXcoord              = 3
    idColumnYcoord              = 4
    idColumnZcoord              = 5
    idColumnThermalGradient     = 15
    idColumnCoolingRate         = 16
    idColumnIdLayer             = 17
    idColumnIdRegion            = 18
    idColumnNumberOfRemelts     = 19
    
    if gradingDistribution == "625centerTo738":
        listFreezingRangeEquilibriumRegions        =   [66.18,66.18,57.836,118.892,113.398,110.245,118.48]
        listFreezingRangeNonEquilibriumRegions     =   [234.33,234.33,209.62,227.94,234.6,227.74,252.19]
        listDiffusivityLiquidRegions               =   [3E-9,3E-9,3E-9,3E-9,3E-9,3E-9,3E-9] # m^2/s
        listCoefficientGibbsThompsonRegions        =   [1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7] # K-m
        listCoefficientSolutePartitionRegions      =   [0.48,0.48,0.48,0.48,0.48,0.48,0.48]
        listCoefficientHarmonicPerturbationRegions =   [28,28,28,28,28,28,28]    
        
        aspectRatio_1           = 0.75
        borderpad_1             = 0.2
        labelspacing_1          = 0.2
        handlelength_1          = 0.4
        handletextpad_1         = 0.2
        borderaxespad_1         = 0.35
        columnspacing_1         = 0.5
        frameon_1               = True
        ncol_1                  = 4
        textPositionFactors_X1  = [0.50,0.25,0.45]
        textPositionFactors_Y1  = [0.95,0.42,0.30]
        rotationAngle_1         = 28
        
        aspectRatio_2           = 0.99
        borderpad_2             = 0.2
        labelspacing_2          = 0.2
        handlelength_2          = 0.4
        handletextpad_2         = 0.2
        borderaxespad_2         = 0.35
        columnspacing_2         = 0.5
        frameon_2               = True
        ncol_2                  = 3
        textPositionFactors_X2  = [0.20,0.15,0.30]
        textPositionFactors_Y2  = [0.95,0.55,0.40]
        rotationAngle_2         = 30
        
    elif gradingDistribution == "738centerTo625":
        listFreezingRangeEquilibriumRegions        =   [100,100,110.245,113.398,118.892,57.836,66.18]
        listFreezingRangeNonEquilibriumRegions     =   [390,390,227.74,234.6,227.94,209.62,234.33]
        listDiffusivityLiquidRegions               =   [3E-9,3E-9,3E-9,3E-9,3E-9,3E-9,3E-9] # m^2/s
        listCoefficientGibbsThompsonRegions        =   [1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7] # K-m
        listCoefficientSolutePartitionRegions      =   [0.48,0.48,0.48,0.48,0.48,0.48,0.48]
        listCoefficientHarmonicPerturbationRegions =   [28,28,28,28,28,28,28]
        
        aspectRatio_1           = 0.75
        borderpad_1             = 0.2
        labelspacing_1          = 0.2
        handlelength_1          = 0.4
        handletextpad_1         = 0.2
        borderaxespad_1         = 0.35
        columnspacing_1         = 0.5
        frameon_1               = True
        ncol_1                  = 4
        textPositionFactors_X1  = [0.50,0.15,0.40]
        textPositionFactors_Y1  = [0.99,0.37,0.30]
        rotationAngle_1         = 23
        
        aspectRatio_2           = 0.99
        borderpad_2             = 0.2
        labelspacing_2          = 0.2
        handlelength_2          = 0.4
        handletextpad_2         = 0.2
        borderaxespad_2         = 0.35
        columnspacing_2         = 0.5
        frameon_2               = True
        ncol_2                  = 3
        textPositionFactors_X2  = [0.20,0.15,0.30]
        textPositionFactors_Y2  = [0.95,0.55,0.40]
        rotationAngle_2         = 30

    # --------- include all the layers
    #listIdLayersToExclude       = []
    #listIdRegionsToExclude      = []
    #plotLayersSeparately        = False
    #plotRegionsSeparately       = True
    #plotRadialDataSeparately    = False
    #tagFigure                   = 'allLayers'
    
    # --------- including only layers 74 to 80
    #listIdLayersToExclude      = 7*[ii for ii in range(1,74)]
    #listIdRegionsToExclude     = 73*[1] + 73*[2] + 73*[3] + 73*[4] + 73*[5] + 73*[6] + 73*[7]
    #plotLayersSeparately       = True
    #plotRegionsSeparately      = False
    #plotRadialDataSeparately   = False    
    #tagFigure                  = 'layer_74_to_80'    
    
    # --------- including only layers 1 to 20
    listIdLayersToExclude      =  7*[ii for ii in range(21,81)]
    listIdRegionsToExclude     = [1 for ii in range(21,81)]  + [2 for ii in range(21,81)]  + [3 for ii in range(21,81)]  + [4 for ii in range(21,81)]  + [5 for ii in range(21,81)]  + [6 for ii in range(21,81)]  + [7 for ii in range(21,81)] 
    plotRegionsSeparately      = False
    plotRadialDataSeparately   = True
    tagFigure                  = 'layer_1_to_20'

    # --------- including only layer 21
    #listIdLayersToExclude      =  7*[ii for ii in range(1,21)] + 7*[ii for ii in range(22,81)]
    #listIdRegionsToExclude     = [1 for ii in range(1,21)]  + [2 for ii in range(1,21)]  + [3 for ii in range(1,21)]  + [4 for ii in range(1,21)]  + [5 for ii in range(1,21)]  + [6 for ii in range(1,21)]  + [7 for ii in range(1,21)] + [1 for ii in range(22,81)]  + [2 for ii in range(22,81)]  + [3 for ii in range(22,81)]  + [4 for ii in range(22,81)]  + [5 for ii in range(22,81)]  + [6 for ii in range(22,81)]  + [7 for ii in range(22,81)] 
    #plotRegionsSeparately      = False
    #plotRadialDataSeparately   = True
    #tagFigure                  = 'layer_21'

    # --------- including only layer 80
    #listIdLayersToExclude      = 7*[ii for ii in range(1,80)]
    #listIdRegionsToExclude     = 79*[1] + 79*[2] + 79*[3] + 79*[4] + 79*[5] + 79*[6] + 79*[7]
    #plotRegionsSeparately      = True
    #plotRadialDataSeparately   = True
    #tagFigure                  = 'layer_80'
    
    # --------- including only layer 79
    #listIdLayersToExclude      = 7*[ii for ii in range(1,79)] + 7*[80]
    #listIdRegionsToExclude     = 78*[1] + 78*[2] + 78*[3] + 78*[4] + 78*[5] + 78*[6] + 78*[7] + [1,2,3,4,5,6,7]
    #plotRegionsSeparately      = True
    #tagFigure                  = 'layer_79'

    # --------- including only layer 78
    #listIdLayersToExclude      = 7*[ii for ii in range(1,78)] + 7*[79] + 7*[80]
    #listIdRegionsToExclude     = 77*[1] + 77*[2] + 77*[3] + 77*[4] + 77*[5] + 77*[6] + 77*[7] + 2*[1,2,3,4,5,6,7]
    #plotRegionsSeparately      = True
    #plotRadialDataSeparately   = True
    #tagFigure                  = 'layer_78'
    
    # --------- including only layer 77
    #listIdLayersToExclude      = 7*[ii for ii in range(1,77)] + 7*[78] + 7*[79] + 7*[80]
    #listIdRegionsToExclude     = 76*[1] + 76*[2] + 76*[3] + 76*[4] + 76*[5] + 76*[6] + 76*[7] + 3*[1,2,3,4,5,6,7]
    #plotRegionsSeparately      = True
    #tagFigure                  = 'layer_77'
    
    # --------- including only layer 76
    #listIdLayersToExclude      = 7*[ii for ii in range(1,76)] + 7*[77] + 7*[78] + 7*[79] + 7*[80]
    #listIdRegionsToExclude     = 75*[1] + 75*[2] + 75*[3] + 75*[4] + 75*[5] + 75*[6] + 75*[7] + 4*[1,2,3,4,5,6,7]
    #plotRegionsSeparately      = True
    #tagFigure                  = 'layer_76'
    
    # --------- including only layer 75
    #listIdLayersToExclude      = 7*[ii for ii in range(1,75)] + 7*[76] + 7*[77] + 7*[78] + 7*[79] + 7*[80]
    #listIdRegionsToExclude     = 74*[1] + 74*[2] + 74*[3] + 74*[4] + 74*[5] + 74*[6] + 74*[7] + 5*[1,2,3,4,5,6,7]
    #plotRegionsSeparately      = True
    #tagFigure                  = 'layer_75'
    
    # --------- including only layer 74
    #listIdLayersToExclude      = 7*[ii for ii in range(1,74)] + 7*[75] + 7*[76] + 7*[77] + 7*[78] + 7*[79] + 7*[80]
    #listIdRegionsToExclude     = 73*[1] + 73*[2] + 73*[3] + 73*[4] + 73*[5] + 73*[6] + 73*[7] + 6*[1,2,3,4,5,6,7]
    #plotRegionsSeparately      = True
    #plotRadialDataSeparately   = True
    #tagFigure                  = 'layer_74'

    # --------- only layers 21 to 60
    #listIdLayersToExclude      = 7*[ii for ii in range(1,21)] + 7*[61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80]
    #listIdRegionsToExclude     = 20*[1] + 20*[2] + 20*[3] + 20*[4] + 20*[5] + 20*[6] + 20*[7] +  20*[1] + 20*[2] + 20*[3] + 20*[4] + 20*[5] + 20*[6] + 20*[7] 
    #plotRegionsSeparately      = True
    #plotRadialDataSeparately   = True
    #tagFigure                  = 'layer_21_to_60'

    # --------- including only layers 1 to 5
    #listIdLayersToExclude      =  7*[ii for ii in range(6,81)]
    #listIdRegionsToExclude     = [1 for ii in range(6,81)]  + [2 for ii in range(6,81)]  + [3 for ii in range(6,81)]  + [4 for ii in range(6,81)]  + [5 for ii in range(6,81)]  + [6 for ii in range(6,81)]  + [7 for ii in range(6,81)] 
    #plotRegionsSeparately      = True
    #plotLayersSeparately       = True
    #plotRadialDataSeparately   = True
    #tagFigure                  = 'layer_1_to_5'
    
    # --------- including only layers 1 to 10
    #listIdLayersToExclude      =  7*[ii for ii in range(11,81)]
    #listIdRegionsToExclude     = [1 for ii in range(11,81)]  + [2 for ii in range(11,81)]  + [3 for ii in range(11,81)]  + [4 for ii in range(11,81)]  + [5 for ii in range(11,81)]  + [6 for ii in range(11,81)]  + [7 for ii in range(11,81)] 
    #plotRegionsSeparately      = True
    #plotLayersSeparately       = True
    #plotRadialDataSeparately   = True
    #tagFigure                  = 'layer_1_to_10'
    
    # --------- including only the regions 1 to 3 of layer 50
    #listIdLayersToExclude      = 7*[ii for ii in range(1,50)] + 7*[ii for ii in range(51,81)] + [50,50,50,50]
    #listIdRegionsToExclude     = 49*[1] + 49*[2] + 49*[3] + 49*[4] + 49*[5] + 49*[6] + 49*[7] + 30*[1] + 30*[2] + 30*[3] + 30*[4] + 30*[5] + 30*[6] + 30*[7] + [4,5,6,7]
    #plotRegionsSeparately      = True
    #tagFigure                  = 'layer_50_regions_1_to_3'
    
    if len(listIdLayersToExclude) != len(listIdRegionsToExclude):
        raise RuntimeError('Number of layers and regions to be excluded must be the same')
    
    listCoordinatesX,listCoordinatesY,listCoordinatesZ,listCoordinatesRadial,listThermalGradients,listCoolingRates,listSolidificationVelocities,listLayerIDs,listRegionsIDs,listNumberOfRemelts=[[] for ii in range(10)]
    listPDASTrivediModel,listPDASKurzFisherModel = ([] for ii in range(2))
    with open(nameFile) as fileCurrent:    
        dataTemp = fileCurrent.readlines()
        
        for iData in range(1, len(dataTemp)):
            dataItems = dataTemp[iData].strip().split()
            coordinateX = float(dataItems[idColumnXcoord-1].strip())
            coordinateY = float(dataItems[idColumnYcoord-1].strip())
            coordinateZ = float(dataItems[idColumnZcoord-1].strip())
            thermalGradient = float(dataItems[idColumnThermalGradient-1].strip())*multiplierThermalGradient
            coolingRate = float(dataItems[idColumnCoolingRate-1].strip())
            solidificationVelocity = coolingRate/thermalGradient
            layerID = int(dataItems[idColumnIdLayer-1].strip())
            regionID = int(dataItems[idColumnIdRegion-1].strip())
            nRemelts = int(dataItems[idColumnNumberOfRemelts-1].strip())
            
            shouldConsiderThisData = True
            
            if listIdLayersToExclude != [] and listIdRegionsToExclude != []:
                for ii in range(0,len(listIdLayersToExclude)):
                    layerIDtoExclude = listIdLayersToExclude[ii]
                    regionIDtoExclude = listIdRegionsToExclude[ii]
                    if layerIDtoExclude == layerID and regionIDtoExclude == regionID:
                        shouldConsiderThisData = False
                        break

            if shouldConsiderThisData == True:
               listCoordinatesX.append(coordinateX)
               listCoordinatesY.append(coordinateY)
               listCoordinatesZ.append(coordinateZ)
               listCoordinatesRadial.append(math.sqrt(coordinateX**2 + coordinateY**2))               
               listThermalGradients.append(thermalGradient)
               listCoolingRates.append(coolingRate)
               listSolidificationVelocities.append(solidificationVelocity)
               listLayerIDs.append(layerID)
               listRegionsIDs.append(regionID)
               listNumberOfRemelts.append(nRemelts)
               coefficientSolutePartition = listCoefficientSolutePartitionRegions[regionID-1]
               freezingRangeEquilibrium = listFreezingRangeEquilibriumRegions[regionID-1]
               freezingRangeNonEquilibrium = listFreezingRangeNonEquilibriumRegions[regionID-1]
               coefficientHarmonicPerturbation = listCoefficientHarmonicPerturbationRegions[regionID-1]
               diffusivityLiquid = listDiffusivityLiquidRegions[regionID-1]
               coefficientGibbsThompson = listCoefficientGibbsThompsonRegions[regionID-1]
               PDAS = computePDAS_NiSuperAlloys(coefficientSolutePartition,freezingRangeEquilibrium,
                                                freezingRangeNonEquilibrium,coefficientHarmonicPerturbation,diffusivityLiquid,
                                                coefficientGibbsThompson,solidificationVelocity,thermalGradient,"Trivedi")
               listPDASTrivediModel.append(PDAS)               
               PDAS = computePDAS_NiSuperAlloys(coefficientSolutePartition,freezingRangeEquilibrium,
                                                freezingRangeNonEquilibrium,coefficientHarmonicPerturbation,diffusivityLiquid,
                                                coefficientGibbsThompson,solidificationVelocity,thermalGradient,"KurzFisher")
               listPDASKurzFisherModel.append(PDAS)
                    
    fileCurrent.close()
    
    # -------------------------------------------------------------------------------------------------------------
    # ------------------------------------- Plot scatter plots in 3D
    # -------------------------------------------------------------------------------------------------------------
    
    if nPointsToPlot > len(listCoordinatesX):
        nPointsToPlot = len(listCoordinatesX)
    
    if listCoordinatesX != []:
        listIDsToPlot = []
        for ii in range(0,nPointsToPlot):
            listIDsToPlot.append(rand.randint(0,len(listCoordinatesX)-1))
        
        listCoordinatesX_toPlot = [listCoordinatesX[ii] for ii in listIDsToPlot]
        listCoordinatesY_toPlot = [listCoordinatesY[ii] for ii in listIDsToPlot]
        listCoordinatesZ_toPlot = [listCoordinatesZ[ii] for ii in listIDsToPlot]
        listThermalGradients_toPlot = [listThermalGradients[ii] for ii in listIDsToPlot]
        listSolidificationVelocities_toPlot = [listSolidificationVelocities[ii] for ii in listIDsToPlot]
        listCoolingRates_toPlot = [listCoolingRates[ii] for ii in listIDsToPlot]

        if plotScatterPlotsIn3D == True:
            plotScatter4D(listCoordinatesX_toPlot,
                          listCoordinatesY_toPlot,
                          listCoordinatesZ_toPlot,
                          listThermalGradients_toPlot, 
                          "Thermal Gradient (K/m)", "X", "Y", "Z")
                          
            plotScatter4D(listCoordinatesX_toPlot, 
                          listCoordinatesY_toPlot, 
                          listCoordinatesZ_toPlot, 
                          listSolidificationVelocities_toPlot, 
                          "Solidification Velocity (m/s)", "X", "Y", "Z")
            
            plotScatter4D(listCoordinatesX_toPlot, 
                          listCoordinatesY_toPlot, 
                          listCoordinatesZ_toPlot, 
                          listCoolingRates_toPlot, 
                          "coolingRates (K/s)", "X", "Y", "Z")     

    # -------------------------------------------------------------------------------------------------------------
    # ------------------------------------- Plot scatter plots in 2D for R and G (Layerwise)
    # -------------------------------------------------------------------------------------------------------------
    
        # Experimental data
        listColumnarR = [0.04732,2.86E-05,2.81E-05,2.48E-05,2.48E-05,2.43E-05,2.39E-05,2.51E-05,2.43E-05]
        listColumnarG = [75151.22313,528.64597,341.75622,109.06193,17.16599,1.93789,1.1529,11.09613,1.22697]
        
        listEquiaxedR = [1.01283,2.81E-05,2.51E-05,2.43E-05]
        listEquiaxedG = [80956.3094,44.61794,11.09613,1.22697]
                
        listColors = ["red","green","blue","black","cyan","magenta","orange","tab:red","tab:green","tab:blue","tab:black","tab:cyan","tab:magenta","tab:orange"]
        listMarkers = ["o","o","o","o","o","o","o","o","o","o","o","o","o","o"]
        #listMarkers = ["o","*","x","s","D","v","^"]       
        if plotRegionsSeparately == True:
            fig = plt.figure()
            ax = fig.add_subplot()
            for iRegion in range(1,8):
                listThermalGradients_ = []
                listSolidificationVelocities_ = []
                for iData in range(0,len(listThermalGradients)):
                    if listRegionsIDs[iData] == iRegion:
                        listThermalGradients_.append(listThermalGradients[iData])
                        listSolidificationVelocities_.append(listSolidificationVelocities[iData])
                        
                plt.scatter(listSolidificationVelocities_,listThermalGradients_, 
                            c = listColors[iRegion-1],label = 'R' + str(iRegion), 
                            alpha = 1, s = 25, marker = listMarkers[iRegion-1])     

            plt.plot(listColumnarR,listColumnarG,c = 'black',linewidth=lineWidth)
            plt.plot(listEquiaxedR,listEquiaxedG,c = 'black',linewidth=lineWidth)
            plt.xscale("log")
            plt.yscale("log")
            xleft, xright = ax.get_xlim()
            ybottom, ytop = ax.get_ylim()            
            ax.set_aspect(aspectRatio_1)            
            #ax.set_aspect('equal', adjustable = 'box')
            leg = ax.legend(loc = 'lower right',prop=legendFont,borderpad=borderpad_1,labelspacing=labelspacing_1,handlelength=handlelength_1,handletextpad=handletextpad_1,borderaxespad=borderaxespad_1,columnspacing=columnspacing_1,frameon=frameon_1,ncol=ncol_1)
            plt.xlabel("R (m/s)",**fontStyleForLabels)
            plt.ylabel("G (K/m)",**fontStyleForLabels)
            
            plotRange_minX = math.log10(min(min(listEquiaxedR),min(listColumnarR),min(listSolidificationVelocities_)))
            plotRange_maxX = math.log10(max(max(listEquiaxedR),max(listColumnarR),max(listSolidificationVelocities_)))
            plotRange_X = plotRange_maxX - plotRange_minX
            
            plotRange_minY = math.log10(min(min(listEquiaxedG),min(listColumnarG),min(listThermalGradients_)))
            plotRange_maxY = math.log10(max(max(listEquiaxedG),max(listColumnarG),max(listThermalGradients_)))
            plotRange_Y = plotRange_maxY - plotRange_minY
            
            plt.text(10**(plotRange_minX + textPositionFactors_X1[0]*plotRange_X),10**(plotRange_minY + textPositionFactors_Y1[0]*plotRange_Y),'Columnar',fontdict=fontDictPlotText)
            plt.text(10**(plotRange_minX + textPositionFactors_X1[1]*plotRange_X),10**(plotRange_minY + textPositionFactors_Y1[1]*plotRange_Y),'Mixed',rotation=rotationAngle_1,fontdict=fontDictPlotText)
            plt.text(10**(plotRange_minX + textPositionFactors_X1[2]*plotRange_X),10**(plotRange_minY + textPositionFactors_Y1[2]*plotRange_Y),'Equiaxed',fontdict=fontDictPlotText)
            
            plt.xticks(fontsize = fontSizeXticks, weight=fontWeightXticks)
            plt.yticks(fontsize = fontSizeYticks, weight=fontWeightYticks)
            for tick in ax.get_xticklabels():
                tick.set_fontname(fontName)
            for tick in ax.get_yticklabels():
                tick.set_fontname(fontName)             
            plt.savefig(pwd + "\\" + "plot_G_R_regionWise_" + str(tagFigure) + imgFormat,bbox_inches='tight',dpi=nPixelsPerInch)

    # -------------------------------------------------------------------------------------------------------------
    # ------------------------------------- Plot scatter plots in 2D for R and G (Regionwise)
    # -------------------------------------------------------------------------------------------------------------
    
        if plotLayersSeparately == True:
            nLayersToPlot = 0
            fig = plt.figure()
            ax = fig.add_subplot()
            for iLayer in range(1,81):
                listThermalGradients_ = []
                listSolidificationVelocities_ = []
                for iData in range(0,len(listThermalGradients)):
                    if listLayerIDs[iData] == iLayer:
                        listThermalGradients_.append(listThermalGradients[iData])
                        listSolidificationVelocities_.append(listSolidificationVelocities[iData])
                
                if listThermalGradients_ != []:
                    nLayersToPlot = nLayersToPlot + 1
                    plt.scatter(listSolidificationVelocities_,listThermalGradients_, 
                                c = listColors[nLayersToPlot - 1], label = 'L' + str(iLayer), 
                                alpha = 1, s = 25, marker = listMarkers[nLayersToPlot - 1])                                                        

            plt.plot(listColumnarR,listColumnarG,c = 'black',linewidth=lineWidth)
            plt.plot(listEquiaxedR,listEquiaxedG,c = 'black',linewidth=lineWidth)
            plt.xscale("log")
            plt.yscale("log")  
            ax.set_aspect(aspectRatio_2)             
            #ax.set_aspect('equal', adjustable = 'box')
            leg = ax.legend(loc = 'lower right',prop=legendFont,borderpad=borderpad_2,labelspacing=labelspacing_2,handlelength=handlelength_2,handletextpad=handletextpad_2,borderaxespad=borderaxespad_2,columnspacing=columnspacing_2,frameon=frameon_2,ncol=ncol_2)
            plt.xlabel("R (m/s)",**fontStyleForLabels)
            plt.ylabel("G (K/m)",**fontStyleForLabels)
            
            plotRange_minX = math.log10(min(min(listEquiaxedR),min(listColumnarR),min(listSolidificationVelocities_)))
            plotRange_maxX = math.log10(max(max(listEquiaxedR),max(listColumnarR),max(listSolidificationVelocities_)))
            plotRange_X = plotRange_maxX - plotRange_minX
            
            plotRange_minY = math.log10(min(min(listEquiaxedG),min(listColumnarG),min(listThermalGradients_)))
            plotRange_maxY = math.log10(max(max(listEquiaxedG),max(listColumnarG),max(listThermalGradients_)))
            plotRange_Y = plotRange_maxY - plotRange_minY
            
            plt.text(10**(plotRange_minX + textPositionFactors_X2[0]*plotRange_X),10**(plotRange_minY + textPositionFactors_Y2[0]*plotRange_Y),'Columnar',fontdict=fontDictPlotText)
            plt.text(10**(plotRange_minX + textPositionFactors_X2[1]*plotRange_X),10**(plotRange_minY + textPositionFactors_Y2[1]*plotRange_Y),'Mixed',rotation=rotationAngle_2,fontdict=fontDictPlotText)
            plt.text(10**(plotRange_minX + textPositionFactors_X2[2]*plotRange_X),10**(plotRange_minY + textPositionFactors_Y2[2]*plotRange_Y),'Equiaxed',fontdict=fontDictPlotText)
            
            plt.xticks(fontsize = fontSizeXticks, weight=fontWeightXticks)
            plt.yticks(fontsize = fontSizeYticks, weight=fontWeightYticks)
            for tick in ax.get_xticklabels():
                tick.set_fontname(fontName)
            for tick in ax.get_yticklabels():
                tick.set_fontname(fontName)            
            plt.savefig(pwd + "\\" + "plot_G_R_layerWise_" + str(tagFigure) + imgFormat,bbox_inches='tight',dpi=nPixelsPerInch)
    
    # ----------------------------------------------------------------------------------------------------------------------
    # ------------------------------------- Plot data as a function of radial distance from the center of the deposit
    # ----------------------------------------------------------------------------------------------------------------------
    
        if plotRadialDataSeparately:
            listCoordinatesRadial_,listThermalGradients_,listSolidificationVelocities_,listCoolingRates_,listPDASTrivediModel_,listPDASKurzFisherModel_ = ([] for ii in range(6))
            listMeanCoordinatesRadial_,listMeanThermalGradients_,listMeanSolidificationVelocities_,listMeanCoolingRates_,listMeanPDASTrivediModel_,listMeanPDASKurzFisherModel_ = ([] for ii in range(6))
            listStdCoordinatesRadial_,listStdThermalGradients_,listStdSolidificationVelocities_,listStdCoolingRates_,listStdPDASTrivediModel_,listStdPDASKurzFisherModel_ = ([] for ii in range(6))

            for iRegion in range(0,7):
                listCoordinatesRadial_.append([])
                listThermalGradients_.append([])
                listCoolingRates_.append([])
                listSolidificationVelocities_.append([])
                listPDASTrivediModel_.append([])
                listPDASKurzFisherModel_.append([])
                for iData in range(0,len(listThermalGradients)):
                    if listRegionsIDs[iData] == iRegion + 1:
                        listCoordinatesRadial_[iRegion].append(listCoordinatesRadial[iData])
                        listThermalGradients_[iRegion].append(listThermalGradients[iData])
                        listSolidificationVelocities_[iRegion].append(listSolidificationVelocities[iData])
                        listCoolingRates_[iRegion].append(listCoolingRates[iData])
                        listPDASTrivediModel_[iRegion].append(listPDASTrivediModel[iData])
                        listPDASKurzFisherModel_[iRegion].append(listPDASKurzFisherModel[iData])
                
                listMeanCoordinatesRadial_.append(stats.mean(listCoordinatesRadial_[iRegion]))
                listStdCoordinatesRadial_.append(stats.stdev(listCoordinatesRadial_[iRegion]))
                listMeanThermalGradients_.append(stats.mean(listThermalGradients_[iRegion]))
                listStdThermalGradients_.append(stats.stdev(listThermalGradients_[iRegion]))
                listMeanSolidificationVelocities_.append(stats.mean(listSolidificationVelocities_[iRegion]))
                listStdSolidificationVelocities_.append(stats.stdev(listSolidificationVelocities_[iRegion]))
                listMeanCoolingRates_.append(stats.mean(listCoolingRates_[iRegion]))
                listStdCoolingRates_.append(stats.stdev(listCoolingRates_[iRegion]))
                listMeanPDASTrivediModel_.append(stats.mean(listPDASTrivediModel_[iRegion]))
                listStdPDASTrivediModel_.append(stats.stdev(listPDASTrivediModel_[iRegion]))
                listMeanPDASKurzFisherModel_.append(stats.mean(listPDASKurzFisherModel_[iRegion]))
                listStdPDASKurzFisherModel_.append(stats.stdev(listPDASKurzFisherModel_[iRegion]))
                
            # --------------- Thermal gradient  
            plotScatter2D(listMeanCoordinatesRadial_,listMeanThermalGradients_,True,pwd + "\\" + "plot_radial_meanG_" + str(tagFigure) + imgFormat,labelX="Radial distance - r (mm)",labelY="Thermal gradient - G (K/m)",nPixelsPerInch=nPixelsPerInch)
            plotScatter2D(listMeanCoordinatesRadial_,listStdThermalGradients_,True,pwd + "\\" + "plot_radial_stdG_" + str(tagFigure) + imgFormat,labelX="Radial distance - r (mm)",labelY="Thermal gradient - G (K/m)",nPixelsPerInch=nPixelsPerInch)
            plotScatter2DwithErrorBars(listMeanCoordinatesRadial_,listMeanThermalGradients_,listStdThermalGradients_,True,pwd + "\\" + "plot_radial_meanstdG_" + str(tagFigure) + imgFormat,labelX="Radial distance - r (mm)",labelY="Thermal gradient - G (K/m)",nPixelsPerInch=nPixelsPerInch)

            # ----------------- Solidification velocity
            plotScatter2D(listMeanCoordinatesRadial_,listMeanSolidificationVelocities_,True,pwd + "\\" + "plot_radial_meanR_" + str(tagFigure) + imgFormat,labelX="Radial distance - (mm)",labelY="Solidification Velocity - R (m/s)",nPixelsPerInch=nPixelsPerInch)
            plotScatter2D(listMeanCoordinatesRadial_,listStdSolidificationVelocities_,True,pwd + "\\" + "plot_radial_stdR_" + str(tagFigure) + imgFormat,labelX="Radial distance - (mm)",labelY="Solidification Velocity - R (m/s)",nPixelsPerInch=nPixelsPerInch)
            plotScatter2DwithErrorBars(listMeanCoordinatesRadial_,listMeanSolidificationVelocities_,listStdSolidificationVelocities_,True,pwd + "\\" + "plot_radial_meanstdR_" + str(tagFigure) + imgFormat,labelX="Radial distance - r (mm)",labelY="Solidification Velocity - R (m/s)",nPixelsPerInch=nPixelsPerInch)
                     
            # ---------------- Cooling rate
            plotScatter2D(listMeanCoordinatesRadial_,listMeanCoolingRates_,True,pwd + "\\" + "plot_radial_meanCoolingRate_" + str(tagFigure) + imgFormat,labelX="Radial distance - (mm)",labelY="Cooling Rate - (K/s)",nPixelsPerInch=nPixelsPerInch)
            plotScatter2D(listMeanCoordinatesRadial_,listStdCoolingRates_,True,pwd + "\\" + "plot_radial_stdCoolingRate_" + str(tagFigure) + imgFormat,labelX="Radial distance - (mm)",labelY="Cooling Rate - (K/s)",nPixelsPerInch=nPixelsPerInch)
            plotScatter2DwithErrorBars(listMeanCoordinatesRadial_,listMeanCoolingRates_,listStdCoolingRates_,True,pwd + "\\" + "plot_radial_meanstdCoolingRate_" + str(tagFigure)+ imgFormat,labelX="Radial distance - (mm)",labelY="Cooling Rate - (K/s)",nPixelsPerInch=nPixelsPerInch)
            
            # ------------------ PDAS (Trivedi model)
            plotScatter2D(listMeanCoordinatesRadial_,listMeanPDASTrivediModel_,True,pwd + "\\" + "plot_radial_meanPDASTrivedi_" + str(tagFigure) + imgFormat,labelX="Radial distance - (mm)",labelY="PDAS - (m)",nPixelsPerInch=nPixelsPerInch)
            plotScatter2D(listMeanCoordinatesRadial_,listStdPDASTrivediModel_,True,pwd + "\\" + "plot_radial_stdPDASTrivedi_" + str(tagFigure) + imgFormat,labelX="Radial distance - (mm)",labelY="PDAS - (m)",nPixelsPerInch=nPixelsPerInch)
            plotScatter2DwithErrorBars(listMeanCoordinatesRadial_,listMeanPDASTrivediModel_,listStdPDASTrivediModel_,True,pwd + "\\" + "plot_radial_meanstdPDASTrivedi_" + str(tagFigure) + imgFormat,labelX="Radial distance - (mm)",labelY="PDAS - (m)",nPixelsPerInch=nPixelsPerInch)

            # ----------------- PDAS (Kurz-Fisher model)            
            plotScatter2D(listMeanCoordinatesRadial_,listMeanPDASKurzFisherModel_,True,pwd + "\\" + "plot_radial_meanPDASKurzFisher_" + str(tagFigure) + imgFormat,labelX="Radial distance - (mm)",labelY="PDAS - (m)",nPixelsPerInch=nPixelsPerInch)
            plotScatter2D(listMeanCoordinatesRadial_,listStdPDASKurzFisherModel_,True,pwd + "\\" + "plot_radial_stdPDASKurzFisher_" + str(tagFigure) + imgFormat,labelX="Radial distance - (mm)",labelY="PDAS - (m)",nPixelsPerInch=nPixelsPerInch)
            plotScatter2DwithErrorBars(listMeanCoordinatesRadial_,listMeanPDASKurzFisherModel_,listStdPDASKurzFisherModel_,True,pwd + "\\" + "plot_radial_meanstdPDASTrivedi_" + str(tagFigure) + imgFormat,labelX="Radial distance - (mm)",labelY="PDAS - (m)",nPixelsPerInch=nPixelsPerInch)


if __name__ == "__main__":
    main()