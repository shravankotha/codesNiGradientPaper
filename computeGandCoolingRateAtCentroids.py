''' This code compuates the therma gradient and cooling rate at the centroids of all the elements for all the time steps 
    given the nodal temperatures at all the time steps. Nodal temperatures corresponding to each abaqus step are extracted from 
    abaqus odb file before runnning this code. The nodal and elemental numbers are renumbered and these files are required before 
    running this code'''


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
    if nArguments != 12:
        text = """Eleven command line arguments are expected: \
                \n\t(1) Nodal coordinates file (renumbered) name including path \
                \n\t(2) Connectivity file (original) name including path \
                \n\t(3) Connectivity file (renumbered) name including path \
                \n\t(4) directory path containing input nodal temperatures files \
                \n\t(5) output directory path \
                \n\t(6) Reference/initial temperature (in the same units as nodal temperature) \
                \n\t(7) liquidus temperature (in the same units as nodal temperature) \
                \n\t(8) layerNumberStart \
                \n\t(9) layerNumberEnd \
                \n\t(10) trajectoryNumberStart \
                \n\t(11) trajectoryNumberEnd"""
                
        raise RuntimeError(text)
    
    time_start = MPI.Wtime()
    nameFileCoordinatesNodal = sys.argv[1]
    fileName_elementConnectivity_orig = sys.argv[2]
    fileName_elementConnectivity = sys.argv[3]
    nodalTempInputDirectoryPath = sys.argv[4]
    outputPath = sys.argv[5]
    refTemperature = float(sys.argv[6])
    temperatureLiquidus = float(sys.argv[7])
    layerNumberStart = int(sys.argv[8])
    layerNumberEnd = int(sys.argv[9])
    trajectoryNumberStart = int(sys.argv[10])
    trajectoryNumberEnd = int(sys.argv[11])

    nTrajectoriesBeforeACoolingStep = 7
    listFileNamesNodalTemperatures = []
    for iLayer in range(layerNumberStart,layerNumberEnd+1):
        for iTrajectory in range(trajectoryNumberStart,trajectoryNumberEnd+1):
            fileNamePrintingStep = nodalTempInputDirectoryPath + '/time_temps_Printing_layer_' + str(iLayer) + '_traj_' + str(iTrajectory) + '.dat'
            fileNameCoolingStep = nodalTempInputDirectoryPath + '/time_temps_Cooling_layer_' + str(iLayer) + '_step_1.dat'
            listFileNamesNodalTemperatures.append(fileNamePrintingStep)
            if iTrajectory == nTrajectoriesBeforeACoolingStep and iLayer != layerNumberEnd:
                listFileNamesNodalTemperatures.append(fileNameCoolingStep)
                    
    #pathDir = os.path.splitdrive(nameFileCoordinatesNodal)[0] + os.path.split(os.path.splitdrive(nameFileCoordinatesNodal)[1])[0] + "/"
    pathDir = outputPath + "/"
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
    
    listCoordinatesCentroidPerProcessor = []
    listCoordinatesNodalElementalPerProcessor = []
    for iElement in range(iElementStart, iElementEnd):        
        listConnectivityElement = listConnectivity[iElement]
        coordinatesCentroid = np.zeros(3)
        listCoordinatesNodalElemental = [[] for ii in range(0,3)]
        for iNode in range(0,8):
            idNode = listConnectivityElement[iNode]
            for iDim in range(0,3):
                listCoordinatesNodalElemental[iDim].append(listCoordinatesNodal[idNode-1][iDim])
                coordinatesCentroid[iDim] = coordinatesCentroid[iDim] + listCoordinatesNodal[idNode-1][iDim]
        coordinatesCentroid = (1/8)*coordinatesCentroid    
        listCoordinatesCentroidPerProcessor.append(coordinatesCentroid)
        listCoordinatesNodalElementalPerProcessor.append(listCoordinatesNodalElemental)
    
    listTimesSolidification,listCoordinatesCentroidal,listCoolingRates,listTemperaturesNodalAtSolidification,listTemperaturesCentroidalAtSolidification,listThermalGradients,listNodesSurroundingThisCentroid,listFileNames = ([] for ii in range(8))
    
    for iFile in range(0,len(listFileNamesNodalTemperatures)):
        start_iFile = MPI.Wtime()
        fileNameNodalTemperatures = listFileNamesNodalTemperatures[iFile]
        start_read = MPI.Wtime()
        if rank == 0: print('processing file : ',fileNameNodalTemperatures, ' from processor : ',rank, flush=True)
        end_read = MPI.Wtime()
        listTimes, listTemperatures = readTemperatures(fileNameNodalTemperatures)
        nElementsCount = 0
        for iElement in range(iElementStart, iElementEnd):
            listConnectivityElement = listConnectivity[iElement]
            listConnectivityElement_orig = listConnectivity_orig[iElement]
            nElementsCount = nElementsCount + 1
            coordinatesCentroid = listCoordinatesCentroidPerProcessor[nElementsCount-1]
            listCoordinatesNodalElemental = listCoordinatesNodalElementalPerProcessor[nElementsCount-1]
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
                if temperatureCentroidal_cur < temperatureLiquidus and temperatureCentroidal_prev > temperatureLiquidus:
                    listFileNames.append(fileNameNodalTemperatures)
                    timeStartSolidification = time_cur - ((temperatureCentroidal_cur-temperatureLiquidus)/(temperatureCentroidal_cur-temperatureCentroidal_prev))*(time_cur-time_prev)
                    time_ratio = (time_cur-timeStartSolidification)/(time_cur-time_prev)
                    listTemperaturesNodalAtSolidification = [temperaturesElemental_cur[i] - (temperaturesElemental_cur[i]-temperaturesElemental_prev[i])*time_ratio for i in range(0,8)]
                    temperatureCentroidalAtSolidification = temperatureCentroidal_cur - (temperatureCentroidal_cur - temperatureCentroidal_prev)*time_ratio
                    listTimesSolidification.append(timeStartSolidification)
                    listCoordinatesCentroidal.append(coordinatesCentroid)
                    listCoolingRates.append((temperatureCentroidal_prev-temperatureCentroidal_cur)/(time_cur-time_prev))
                    listTemperaturesCentroidalAtSolidification.append(temperatureCentroidalAtSolidification)
                    listNodesSurroundingThisCentroid.append(listConnectivityElement_orig)
                    dN_dxyz = dShapeFunction_dCaterianCoords(listCoordinatesNodalElemental, [0,0,0])
                    magnitude = 0
                    for iDimension in range(0,3):
                        gradient_ = 0
                        for iNode in range(0,8):
                            gradient_ = gradient_ + dN_dxyz[iNode][iDimension]*listTemperaturesNodalAtSolidification[iNode]
                        magnitude = magnitude + gradient_**2
                    listThermalGradients.append(math.sqrt(magnitude))
                temperatureCentroidal_prev = temperatureCentroidal_cur
                temperaturesElemental_prev = temperaturesElemental_cur
                time_prev = time_cur
        end_iFile = MPI.Wtime()        
        if rank == 0: print('File read time : ', str(end_read-start_read), ' Calculation time for all the elements : ',str(end_iFile-start_iFile), ' Total time so far : ', str(MPI.Wtime()-time_start),flush=True)
    
    comm.barrier()    
    print('finished calculations from processor : ',rank, flush=True)
        
    if rank != 0:
        comm.send(listFileNames,dest=0)
        comm.send(listTimesSolidification,dest=0)
        comm.send(listCoordinatesCentroidal,dest=0)
        comm.send(listNodesSurroundingThisCentroid,dest=0)
        comm.send(listCoolingRates,dest=0)
        comm.send(listTemperaturesCentroidalAtSolidification,dest=0)
        comm.send(listThermalGradients,dest=0)
    
    if rank == 0:
        listFileNamesAll = []
        listFileNamesAll.append(listFileNames)
        if size_comm > 0:
            for iRank in range(1,size_comm):
                data = comm.recv(source=iRank)
                listFileNamesAll.append(data)
    
        listTimesSolidificationAll = []
        listTimesSolidificationAll.append(listTimesSolidification)
        if size_comm > 0:
            for iRank in range(1,size_comm):
                data = comm.recv(source=iRank)
                listTimesSolidificationAll.append(data)
            
        listCoordinatesCentroidalAll = []
        listCoordinatesCentroidalAll.append(listCoordinatesCentroidal)
        if size_comm > 0:
            for iRank in range(1,size_comm):
                data = comm.recv(source=iRank)
                listCoordinatesCentroidalAll.append(data)
                
        listNodesSurroundingThisCentroidAll = []
        listNodesSurroundingThisCentroidAll.append(listNodesSurroundingThisCentroid)
        if size_comm > 0:
            for iRank in range(1,size_comm):
                data = comm.recv(source=iRank)
                listNodesSurroundingThisCentroidAll.append(data)
    
        listCoolingRatesAll = []
        listCoolingRatesAll.append(listCoolingRates)
        if size_comm > 0:
            for iRank in range(1,size_comm):
                data = comm.recv(source=iRank)
                listCoolingRatesAll.append(data)
            
        listTemperaturesCentroidalAtSolidificationAll = []
        listTemperaturesCentroidalAtSolidificationAll.append(listTemperaturesCentroidalAtSolidification)
        if size_comm > 0:
            for iRank in range(1,size_comm):
                data = comm.recv(source=iRank)
                listTemperaturesCentroidalAtSolidificationAll.append(data)

        listThermalGradientsAll = []
        listThermalGradientsAll.append(listThermalGradients)
        if size_comm > 0:
            for iRank in range(1,size_comm):
                data = comm.recv(source=iRank)
                listThermalGradientsAll.append(data) 

        out_path = pathDir + 'G_coolingRates_temperatures.dat'
        with open(out_path, 'w') as file_out:
            file_out.write('fileName    stepTime    centroidcoordX  centroidcoordY   centroidcoordZ    surrNode1    surrNode2   surrNode3   surrNode4   surrNode5   surrNode6   surrNode7   surrNode8   tempCentroid    thermalGradient     coolingRate\n')
            for iData in range(0,len(listThermalGradientsAll)):
                for iSubData in range(0,len(listThermalGradientsAll[iData])):
                    file_out.write("{0:45s}{1:25.10f}{2:25.10f}{3:25.10f}{4:25.10f}{5:10d}{6:10d}{7:10d}{8:10d}{9:10d}{10:10d}{11:10d}{12:10d}{13:25.10f}{14:25.10f}{15:25.10f}\n".format(listFileNamesAll[iData][iSubData],
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
                                                                             listCoolingRatesAll[iData][iSubData]
                                                                            ))
        file_out.close() 
        
    if rank == 0 : print('Total time elapsed: ', str(MPI.Wtime()-time_start), ' (s)')
            
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