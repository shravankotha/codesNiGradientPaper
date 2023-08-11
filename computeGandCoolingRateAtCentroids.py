''' This code computes the therma gradient and cooling rate at the centroids of all the elements for all the time steps 
    given the nodal temperatures at all the time steps. Nodal temperatures corresponding to each abaqus step are extracted from 
    abaqus odb file before runnning this code. The nodal and elemental numbers are renumbered and these files are required before 
    running this code
'''

import sys
import os
from mpi4py import MPI
import statistics as stats
import random as rand
import time
import numpy as np
import math
import inspect
from performLinearPartitioning import performLinearPartitioning
from hexElement import dShapeFunction_dCaterianCoords, getJacobianHexElement
from hexElement import dShapeFunction_dNaturalCoords, coordsNaturalElementHex
from readConnectivityFromAFile import readConnectivity
from readNodalCoordsFromAFile import readCoordinatesNodal

def main():
    comm = MPI.COMM_WORLD
    size_comm = comm.Get_size()
    rank = comm.Get_rank()
    nArguments = len(sys.argv)
    if nArguments != 16:
        text = """15 command line arguments are expected: \
                \n\t(1) Nodal coordinates file (renumbered) \
                \n\t(2) Connectivity file (original) \
                \n\t(3) Connectivity file (renumbered) \
                \n\t(4) Reference/initial temperature (in the same units as nodal temperature) \
                \n\t(5) Flag to determine interface (S: for Solidus, L: for Liquidus, A: for avg of Solidus and Liquidus) \
                \n\t(6) Solidus and liquidus scaling factor used in Abaqus simulations for deposit\
                \n\t(7) Solidus and liquidus scaling factor used in Abaqus simulations for base plate\
                \n\t(8) gradingDistribution (625centerTo738/738centerTo625) \
                \n\t(9) layerNumberStart \
                \n\t(10) layerNumberEnd \
                \n\t(11) trajectoryNumberStart \
                \n\t(12) trajectoryNumberEnd \
                \n\t(13) layer thickness (same units as nodal coordinates) \
                \n\t(14) total number of layers in the build \
                \n\t(15) shouldIncludeTheFinalCoolingSteps (0/1)"""
                
        raise RuntimeError(text)
    
    # -----------------------------------------------------------------------------------------------------------------------
    # ------------------ parse all the inputs and options
    # -----------------------------------------------------------------------------------------------------------------------
    
    time_start = MPI.Wtime()
    nameFileCoordinatesNodal = sys.argv[1]
    fileName_elementConnectivity_orig = sys.argv[2]
    outputPath = sys.argv[3]
    refTemperature = float(sys.argv[4])
    interfaceTemperatureFlag = str(sys.argv[5])
    solidusLiquidusScalingFactorDeposit = float(sys.argv[6])
    solidusLiquidusScalingFactorBasePlate = float(sys.argv[7])
    gradingDistribution =  str(sys.argv[8])
    layerNumberStart = int(sys.argv[9])
    layerNumberEnd = int(sys.argv[10])
    trajectoryNumberStart = int(sys.argv[11])
    trajectoryNumberEnd = int(sys.argv[12])
    thicknessLayer = float(sys.argv[13])
    totalNumberOfLayerInTheBuild = int(sys.argv[14])
    shouldIncludeTheFinalCoolingSteps = int(sys.argv[15])
    
    listRadiiRegionsStart = [0,1.713525,2.189868,2.736886,3.22405,3.743497,4.249974]
    listRadiiRegionsEnd = [1.713525,2.189868,2.736886,3.22405,3.743497,4.249974,4.503479]
    if gradingDistribution == "625centerTo738":
        listSolidusTemperaturesRegionsOriginal = [1126,1126,1147,1125,1115,1118,1088] # equilibrium -- [1290,1290,1298,1251,1246,1252,1277]   
        listLiquidusTemperaturesRegionsOriginal = [1360,1360,1357,1353,1349,1345,1340] # equilibrium -- [1350,1350,1379,1346,1348,1353,1337]  
    elif gradingDistribution == "738centerTo625":
        listSolidusTemperaturesRegionsOriginal = [1088,1088,1118,1115,1125,1147,1126]    # equilibrium -- [1277,1277,1252,1246,1251,1298,1290]
        listLiquidusTemperaturesRegionsOriginal = [1340,1340,1345,1349,1353,1357,1360]   # equilibrium -- [1337,1337,1353,1348,1346,1379,1350] 
    
    # -----------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------
    
    scalingFactor = solidusLiquidusScalingFactorDeposit
    scalingFactorBasePlate = solidusLiquidusScalingFactorBasePlate
    # scale solidus and liquidus to be consistent with abaqus 
    listSolidusTemperaturesRegions, listLiquidusTemperaturesRegions = [[] for ii in range(2)]
    for ii in range(0,len(listSolidusTemperaturesRegionsOriginal)):        
        Tsolidus_new = 0.5*(listSolidusTemperaturesRegionsOriginal[ii]+listLiquidusTemperaturesRegionsOriginal[ii] - scalingFactor*(listLiquidusTemperaturesRegionsOriginal[ii]-listSolidusTemperaturesRegionsOriginal[ii]))
        Tliquidus_new = 0.5*(listSolidusTemperaturesRegionsOriginal[ii]+listLiquidusTemperaturesRegionsOriginal[ii] + scalingFactor*(listLiquidusTemperaturesRegionsOriginal[ii]-listSolidusTemperaturesRegionsOriginal[ii]))
        listSolidusTemperaturesRegions.append(Tsolidus_new)
        listLiquidusTemperaturesRegions.append(Tliquidus_new)        
    listAverageTemperaturesRegions = [0.5*(listSolidusTemperaturesRegions[i]+listLiquidusTemperaturesRegions[i]) for i in range(0,len(listLiquidusTemperaturesRegions))]
    solidusTemperatureBasePlateOriginal = 1125      # equilibrium -- 1440
    liquidusTemperatureBasePlateOriginal = 1518     # equilibrium -- 1505
    solidusTemperatureBasePlate = 0.5*(solidusTemperatureBasePlateOriginal+liquidusTemperatureBasePlateOriginal - scalingFactorBasePlate*(liquidusTemperatureBasePlateOriginal-solidusTemperatureBasePlateOriginal))
    liquidusTemperatureBasePlate = 0.5*(solidusTemperatureBasePlateOriginal+liquidusTemperatureBasePlateOriginal + scalingFactorBasePlate*(liquidusTemperatureBasePlateOriginal-solidusTemperatureBasePlateOriginal))
    averageTemperatureBasePlate = 0.5*(solidusTemperatureBasePlate + liquidusTemperatureBasePlate)    

    nTrajectoriesBeforeACoolingStep = 7
    
    #
    listFileNamesNodalTemperatures = []
    for iLayer in range(layerNumberStart,layerNumberEnd+1):
        for iTrajectory in range(trajectoryNumberStart,trajectoryNumberEnd+1):
            fileNamePrintingStep = 'time_temps_Printing_layer_' + str(iLayer) + '_traj_' + str(iTrajectory) + '.dat'
            fileNameCoolingStep = 'time_temps_Cooling_layer_' + str(iLayer) + '_step_1.dat'
            listFileNamesNodalTemperatures.append(fileNamePrintingStep)
            if iTrajectory == nTrajectoriesBeforeACoolingStep and iLayer != layerNumberEnd:
                listFileNamesNodalTemperatures.append(fileNameCoolingStep)
    if shouldIncludeTheFinalCoolingSteps == 1:
       listFileNamesNodalTemperatures.append('time_temps_Cooling.dat')
       listFileNamesNodalTemperatures.append('time_temps_MoreCooling.dat')
       
    #pathDir = os.path.splitdrive(nameFileCoordinatesNodal)[0] + os.path.split(os.path.splitdrive(nameFileCoordinatesNodal)[1])[0] + "/"
    #pathDir = outputPath + "/"
    # -------------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------- Read nodal coordinates and connectivity
    # -------------------------------------------------------------------------------------------------------------------------------
    listCoordinatesNodal = readCoordinatesNodal(nameFileCoordinatesNodal)
    listConnectivity = readConnectivity(fileName_elementConnectivity)
    listConnectivity_orig = readConnectivity(fileName_elementConnectivity_orig)
    # -------------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------- Read connectivity information
    # -------------------------------------------------------------------------------------------------------------------------------
    startIndicesPartition, endIndicesPartition = performLinearPartitioning(len(listConnectivity),size_comm)
    iElementStart = startIndicesPartition[rank]
    iElementEnd = endIndicesPartition[rank]
    nElementsPerProcessor = iElementEnd - iElementStart + 1

    
    if interfaceTemperatureFlag == "S":
        listTemperaturesInterfaceRegions = listSolidusTemperaturesRegions
        temperatureInterfaceBasePlate = solidusTemperatureBasePlate
    elif interfaceTemperatureFlag == "L":
        listTemperaturesInterfaceRegions = listLiquidusTemperaturesRegions
        temperatureInterfaceBasePlate = liquidusTemperatureBasePlate
    elif interfaceTemperatureFlag == "A":
        listTemperaturesInterfaceRegions = listAverageTemperaturesRegions
        temperatureInterfaceBasePlate = averageTemperatureBasePlate
    else:
        raise RuntimeError('Flag to determine the interface is not valid')
            
    listCoordinatesCentroidPerProcessor = []
    listCoordinatesNodalElementalPerProcessor = []
    listTemperaturesInterfaceElementalPerProcessor = []
    listRegionIDsElementalPerProcessor = []
    listLayerIDsElementalPerProcessor = []
    for iElement in range(iElementStart, iElementEnd + 1):
        listConnectivityElement = listConnectivity[iElement]
        coordinatesCentroid = np.zeros(3)
        listCoordinatesNodalElemental = [[] for ii in range(0,3)]
        for iNode in range(0,8):
            idNode = listConnectivityElement[iNode]
            for iDim in range(0,3):
                listCoordinatesNodalElemental[iDim].append(listCoordinatesNodal[idNode-1][iDim])
                coordinatesCentroid[iDim] = coordinatesCentroid[iDim] + listCoordinatesNodal[idNode-1][iDim]
        coordinatesCentroid = (1/8)*coordinatesCentroid
        radiiCentroid = math.sqrt(coordinatesCentroid[0]**2 + coordinatesCentroid[1]**2)
        if coordinatesCentroid[2] > 0:
            for iRegion in range(0,len(listRadiiRegionsStart)):
                if radiiCentroid > listRadiiRegionsStart[iRegion] and radiiCentroid < listRadiiRegionsEnd[iRegion]:
                    listTemperaturesInterfaceElementalPerProcessor.append(listTemperaturesInterfaceRegions[iRegion])
                    listRegionIDsElementalPerProcessor.append(iRegion + 1)
            for layerID in range(0,totalNumberOfLayerInTheBuild):
                if coordinatesCentroid[2] > thicknessLayer*layerID and coordinatesCentroid[2] <= thicknessLayer*(layerID+1):
                    listLayerIDsElementalPerProcessor.append(layerID+1)
        else:
             listTemperaturesInterfaceElementalPerProcessor.append(temperatureInterfaceBasePlate)
             listRegionIDsElementalPerProcessor.append(0)
             listLayerIDsElementalPerProcessor.append(0)
                
        listCoordinatesCentroidPerProcessor.append(coordinatesCentroid)
        listCoordinatesNodalElementalPerProcessor.append(listCoordinatesNodalElemental)
    
    print('rank:',rank,'nElements,startingID,endingID:',len(listTemperaturesInterfaceElementalPerProcessor), startIndicesPartition[rank]+1, endIndicesPartition[rank]+1)
    comm.barrier()
    listTimesSolidification,listCoordinatesCentroidal,listCoolingRates,listTemperaturesNodalAtSolidification, \
    listTemperaturesCentroidalAtSolidification,listThermalGradients,listNodesSurroundingThisCentroid,listFileNames,listRegionIDs, listLayerIDs = ([] for ii in range(10))
    
    listTempPerProcessor = [0 for i in range(nElementsPerProcessor)]
    listTimesSolidificationConsideringRemelting = [0 for i in range(nElementsPerProcessor)]
    listCoordinatesCentroidalConsideringRemelting = [[0,0,0] for i in range(nElementsPerProcessor)]
    listCoolingRatesConsideringRemelting = [0 for i in range(nElementsPerProcessor)]
    listTemperaturesNodalAtSolidificationConsideringRemelting = [0 for i in range(nElementsPerProcessor)]
    listTemperaturesCentroidalAtSolidificationConsideringRemelting = [0 for i in range(nElementsPerProcessor)]
    listThermalGradientsConsideringRemelting = [0 for i in range(nElementsPerProcessor)]
    listNodesSurroundingThisCentroidConsideringRemelting = [[0,0,0,0,0,0,0,0] for i in range(nElementsPerProcessor)]
    listFileNamesConsideringRemelting = [0 for i in range(nElementsPerProcessor)]
    listNumberOfRemelts = [0 for i in range(nElementsPerProcessor)]
    listRegionIDsConsideringRemelting = [-1 for i in range(nElementsPerProcessor)]
    listLayerIDsConsideringRemelting =  [-1 for i in range(nElementsPerProcessor)]
     
    for iFile in range(0,len(listFileNamesNodalTemperatures)):
        start_iFile = MPI.Wtime()
        fileNameNodalTemperatures = listFileNamesNodalTemperatures[iFile]
        start_read = MPI.Wtime()
        if rank == 0: print('processing file : ',fileNameNodalTemperatures, ' from processor : ',rank)
        end_read = MPI.Wtime()
        listTimes, listTemperatures = readTemperatures(fileNameNodalTemperatures)
        nElementsCount = 0
        for iElement in range(iElementStart, iElementEnd+1):
            listConnectivityElement = listConnectivity[iElement]
            listConnectivityElement_orig = listConnectivity_orig[iElement]
            nElementsCount = nElementsCount + 1
            coordinatesCentroid = listCoordinatesCentroidPerProcessor[nElementsCount-1]            
            listCoordinatesNodalElemental = listCoordinatesNodalElementalPerProcessor[nElementsCount-1]
            listTemperatureInterfaceElement = listTemperaturesInterfaceElementalPerProcessor[nElementsCount-1]
            temperaturesElemental_prev = [refTemperature for i in range(0,8)]
            temperatureCentroidal_prev = refTemperature
            time_prev = 0
            for iTime in range(0,len(listTimes)):
                time_cur = float(listTimes[iTime])
                temperaturesElemental_cur = []
                for iNode in range(0,8):
                    idNode = listConnectivityElement[iNode]
                    temperaturesElemental_cur.append(float(listTemperatures[iTime][idNode-1]))
                temperatureCentroidal_cur = stats.mean(temperaturesElemental_cur)                
                if temperatureCentroidal_cur < listTemperatureInterfaceElement and temperatureCentroidal_prev > listTemperatureInterfaceElement:

                    timeStartSolidification = time_cur - ((temperatureCentroidal_cur-listTemperatureInterfaceElement)/(temperatureCentroidal_cur-temperatureCentroidal_prev))*(time_cur-time_prev)
                    time_ratio = (time_cur-timeStartSolidification)/(time_cur-time_prev)
                    listTemperaturesNodalAtSolidification = [temperaturesElemental_cur[i] - (temperaturesElemental_cur[i]-temperaturesElemental_prev[i])*time_ratio for i in range(0,8)]
                    temperatureCentroidalAtSolidification = temperatureCentroidal_cur - (temperatureCentroidal_cur - temperatureCentroidal_prev)*time_ratio
                    coolingRate = (temperatureCentroidal_prev-temperatureCentroidal_cur)/(time_cur-time_prev)
                    dN_dxyz = dShapeFunction_dCaterianCoords(listCoordinatesNodalElemental, [0,0,0])
                    magnitude = 0
                    for iDimension in range(0,3):
                        gradient_ = 0
                        for iNode in range(0,8):
                            gradient_ = gradient_ + dN_dxyz[iNode][iDimension]*listTemperaturesNodalAtSolidification[iNode]
                        magnitude = magnitude + gradient_**2                    
                    magnitude = math.sqrt(magnitude)                    
                    
                    listFileNames.append(fileNameNodalTemperatures)                    
                    listTimesSolidification.append(timeStartSolidification)
                    listCoordinatesCentroidal.append(coordinatesCentroid)
                    listCoolingRates.append(coolingRate)
                    listTemperaturesCentroidalAtSolidification.append(temperatureCentroidalAtSolidification)
                    listNodesSurroundingThisCentroid.append(listConnectivityElement_orig)
                    listThermalGradients.append(magnitude)
                    listRegionIDs.append(listRegionIDsElementalPerProcessor[nElementsCount-1])
                    listLayerIDs.append(listLayerIDsElementalPerProcessor[nElementsCount-1])
                    
                    listLayerIDsConsideringRemelting[nElementsCount-1] = listLayerIDsElementalPerProcessor[nElementsCount-1]
                    listRegionIDsConsideringRemelting[nElementsCount-1] = listRegionIDsElementalPerProcessor[nElementsCount-1]
                    listNumberOfRemelts[nElementsCount-1] = listNumberOfRemelts[nElementsCount-1] + 1 
                    listFileNamesConsideringRemelting[nElementsCount-1] = fileNameNodalTemperatures
                    listTimesSolidificationConsideringRemelting[nElementsCount-1] = timeStartSolidification
                    listCoordinatesCentroidalConsideringRemelting[nElementsCount-1] = coordinatesCentroid
                    listCoolingRatesConsideringRemelting[nElementsCount-1] = coolingRate
                    listTemperaturesCentroidalAtSolidificationConsideringRemelting[nElementsCount-1] = temperatureCentroidalAtSolidification
                    listNodesSurroundingThisCentroidConsideringRemelting[nElementsCount-1] = listConnectivityElement_orig
                    listThermalGradientsConsideringRemelting[nElementsCount-1] = magnitude
                    
                temperatureCentroidal_prev = temperatureCentroidal_cur
                temperaturesElemental_prev = temperaturesElemental_cur
                time_prev = time_cur
        end_iFile = MPI.Wtime()        
        if rank == 0: print('File read time : ', str(end_read-start_read), ' Calculation time for all the elements : ',str(end_iFile-start_iFile), ' Total time so far : ', str(MPI.Wtime()-time_start))
    
    comm.barrier()    
    print('finished calculations from processor : ',rank)
        
    if rank != 0:
        comm.send(listFileNames,dest=0)
        comm.send(listTimesSolidification,dest=0)
        comm.send(listCoordinatesCentroidal,dest=0)
        comm.send(listNodesSurroundingThisCentroid,dest=0)
        comm.send(listCoolingRates,dest=0)
        comm.send(listTemperaturesCentroidalAtSolidification,dest=0)
        comm.send(listThermalGradients,dest=0)
        comm.send(listRegionIDs,dest=0)
        comm.send(listLayerIDs,dest=0)
        
        comm.send(listLayerIDsConsideringRemelting,dest=0)
        comm.send(listRegionIDsConsideringRemelting,dest=0)
        comm.send(listNumberOfRemelts,dest=0)        
        comm.send(listFileNamesConsideringRemelting,dest=0)
        comm.send(listTimesSolidificationConsideringRemelting,dest=0)
        comm.send(listCoordinatesCentroidalConsideringRemelting,dest=0)
        comm.send(listNodesSurroundingThisCentroidConsideringRemelting,dest=0)
        comm.send(listCoolingRatesConsideringRemelting,dest=0)
        comm.send(listTemperaturesCentroidalAtSolidificationConsideringRemelting,dest=0)
        comm.send(listThermalGradientsConsideringRemelting,dest=0)        
    
    if rank == 0:    
        listFileNamesAll = receiveAndCombineDataFromProcessors([listFileNames],comm)
        listTimesSolidificationAll = receiveAndCombineDataFromProcessors([listTimesSolidification],comm)
        listCoordinatesCentroidalAll = receiveAndCombineDataFromProcessors([listCoordinatesCentroidal],comm)
        listNodesSurroundingThisCentroidAll = receiveAndCombineDataFromProcessors([listNodesSurroundingThisCentroid],comm)
        listCoolingRatesAll = receiveAndCombineDataFromProcessors([listCoolingRates],comm)
        listTemperaturesCentroidalAtSolidificationAll = receiveAndCombineDataFromProcessors([listTemperaturesCentroidalAtSolidification],comm)
        listThermalGradientsAll = receiveAndCombineDataFromProcessors([listThermalGradients],comm)
        listRegionIDsAll = receiveAndCombineDataFromProcessors([listRegionIDs],comm)
        listLayerIDsAll = receiveAndCombineDataFromProcessors([listLayerIDs],comm)
        #
        listLayerIDsConsideringRemeltingAll = receiveAndCombineDataFromProcessors([listLayerIDsConsideringRemelting],comm)
        listRegionIDsConsideringRemeltingAll = receiveAndCombineDataFromProcessors([listRegionIDsConsideringRemelting],comm)
        listNumberOfRemeltsAll = receiveAndCombineDataFromProcessors([listNumberOfRemelts],comm)
        listFileNamesConsideringRemeltingAll = receiveAndCombineDataFromProcessors([listFileNamesConsideringRemelting],comm)
        listTimesSolidificationConsideringRemeltingAll = receiveAndCombineDataFromProcessors([listTimesSolidificationConsideringRemelting],comm)
        listCoordinatesCentroidalConsideringRemeltingAll = receiveAndCombineDataFromProcessors([listCoordinatesCentroidalConsideringRemelting],comm)
        listNodesSurroundingThisCentroidConsideringRemeltingAll = receiveAndCombineDataFromProcessors([listNodesSurroundingThisCentroidConsideringRemelting],comm)
        listCoolingRatesConsideringRemeltingAll = receiveAndCombineDataFromProcessors([listCoolingRatesConsideringRemelting],comm)
        listTemperaturesCentroidalAtSolidificationConsideringRemeltingAll = receiveAndCombineDataFromProcessors([listTemperaturesCentroidalAtSolidificationConsideringRemelting],comm)
        listThermalGradientsConsideringRemeltingAll = receiveAndCombineDataFromProcessors([listThermalGradientsConsideringRemelting],comm)

        out_path = 'G_K_noRemelting_layer_' + str(layerNumberStart) + '_' + str(layerNumberEnd) + '.dat'
        with open(out_path, 'w') as file_out:
            file_out.write('fileName    stepTime    centroidcoordX  centroidcoordY   centroidcoordZ    surrNode1    surrNode2   surrNode3   surrNode4   surrNode5   surrNode6   surrNode7   surrNode8   tempCentroid    thermalGradient     coolingRate     layerID(0 for baseplate)   regionID(0 for baseplate)\n')
            for iData in range(0,len(listThermalGradientsAll)):
                for iSubData in range(0,len(listThermalGradientsAll[iData])):
                    file_out.write("{0:90s}{1:25.10f}{2:25.10f}{3:25.10f}{4:25.10f}{5:10d}{6:10d}{7:10d}{8:10d}{9:10d}{10:10d}{11:10d}{12:10d}{13:25.10f}{14:25.10f}{15:25.10f}{16:10d}{17:10d}\n".format(listFileNamesAll[iData][iSubData],
                                                                             listTimesSolidificationAll[iData][iSubData], 
                                                                             listCoordinatesCentroidalAll[iData][iSubData][0],
                                                                             listCoordinatesCentroidalAll[iData][iSubData][1],
                                                                             listCoordinatesCentroidalAll[iData][iSubData][2],
                                                                             listNodesSurroundingThisCentroidAll[iData][iSubData][0],
                                                                             listNodesSurroundingThisCentroidAll[iData][iSubData][1],
                                                                             listNodesSurroundingThisCentroidAll[iData][iSubData][2],
                                                                             listNodesSurroundingThisCentroidAll[iData][iSubData][3],
                                                                             listNodesSurroundingThisCentroidAll[iData][iSubData][4],
                                                                             listNodesSurroundingThisCentroidAll[iData][iSubData][5],
                                                                             listNodesSurroundingThisCentroidAll[iData][iSubData][6],
                                                                             listNodesSurroundingThisCentroidAll[iData][iSubData][7],                                                                             
                                                                             listTemperaturesCentroidalAtSolidificationAll[iData][iSubData],
                                                                             listThermalGradientsAll[iData][iSubData],
                                                                             listCoolingRatesAll[iData][iSubData],
                                                                             listLayerIDsAll[iData][iSubData],
                                                                             listRegionIDsAll[iData][iSubData]
                                                                            ))
        file_out.close() 
        
        
        out_path = 'G_K_withRemelting_layer_' + str(layerNumberStart) + '_' + str(layerNumberEnd) + '.dat'
        with open(out_path, 'w') as file_out:
            file_out.write('fileName    stepTime    centroidcoordX  centroidcoordY   centroidcoordZ    surrNode1    surrNode2   surrNode3   surrNode4   surrNode5   surrNode6   surrNode7   surrNode8   tempCentroid    thermalGradient     coolingRate     layerID(0 for baseplate)    regionID(0 for baseplate)    nRemelts\n')
            for iData in range(0,len(listThermalGradientsConsideringRemeltingAll)):
                for iSubData in range(0,len(listThermalGradientsConsideringRemeltingAll[iData])):
                    if listTemperaturesCentroidalAtSolidificationConsideringRemeltingAll[iData][iSubData] != 0:
                        file_out.write("{0:90s}{1:25.10f}{2:25.10f}{3:25.10f}{4:25.10f}{5:10d}{6:10d}{7:10d}{8:10d}{9:10d}{10:10d}{11:10d}{12:10d}{13:25.10f}{14:25.10f}{15:25.10f}{16:10d}{17:10d}{18:10d}\n".format(listFileNamesConsideringRemeltingAll[iData][iSubData],
                                                                                listTimesSolidificationConsideringRemeltingAll[iData][iSubData], 
                                                                                listCoordinatesCentroidalConsideringRemeltingAll[iData][iSubData][0],
                                                                                listCoordinatesCentroidalConsideringRemeltingAll[iData][iSubData][1],
                                                                                listCoordinatesCentroidalConsideringRemeltingAll[iData][iSubData][2],
                                                                                listNodesSurroundingThisCentroidConsideringRemeltingAll[iData][iSubData][0],
                                                                                listNodesSurroundingThisCentroidConsideringRemeltingAll[iData][iSubData][1],
                                                                                listNodesSurroundingThisCentroidConsideringRemeltingAll[iData][iSubData][2],
                                                                                listNodesSurroundingThisCentroidConsideringRemeltingAll[iData][iSubData][3],
                                                                                listNodesSurroundingThisCentroidConsideringRemeltingAll[iData][iSubData][4],
                                                                                listNodesSurroundingThisCentroidConsideringRemeltingAll[iData][iSubData][5],
                                                                                listNodesSurroundingThisCentroidConsideringRemeltingAll[iData][iSubData][6],
                                                                                listNodesSurroundingThisCentroidConsideringRemeltingAll[iData][iSubData][7],                                                                             
                                                                                listTemperaturesCentroidalAtSolidificationConsideringRemeltingAll[iData][iSubData],
                                                                                listThermalGradientsConsideringRemeltingAll[iData][iSubData],
                                                                                listCoolingRatesConsideringRemeltingAll[iData][iSubData],
                                                                                listLayerIDsConsideringRemeltingAll[iData][iSubData],
                                                                                listRegionIDsConsideringRemeltingAll[iData][iSubData],
                                                                                listNumberOfRemeltsAll[iData][iSubData]
                                                                                ))
        file_out.close()
        
        
    if rank == 0 : print('Total time elapsed: ', str(MPI.Wtime()-time_start), ' (s)')
            

def receiveAndCombineDataFromProcessors(listDataAll,comm):
    size_comm = comm.Get_size()
    if size_comm > 0:
        for iRank in range(1,size_comm):
            data = comm.recv(source=iRank)
            listDataAll.append(data)            
    return listDataAll        
            
def readTemperatures(nameFile):    
    listTimes, listTemperatures = ([] for ii in range(2))
    with open(nameFile) as fileCurrent:
        dataTemp = fileCurrent.readlines()
        iTimeStep = 0
        for iData in range(0, len(dataTemp)):     
            if str(dataTemp[iData].strip())[0:4] == "time":
                listTimes.append(float(dataTemp[iData].strip().split(":")[1].strip()))
                listTemperatures.append([])
                iTimeStep = iTimeStep + 1                
            else:
                listTemperatures[iTimeStep-1].append(float(dataTemp[iData].strip()))    
    return listTimes, listTemperatures
    
if __name__ == "__main__":
    main()