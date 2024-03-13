'''
 This code is written to extract output variables for a given node.
 node ID is both supplied as cmd line arguments
 This code is element agnostic, tested and works. May be modified to include state variables.
 name of the varible to be extracted can be given as an argument.
'''
from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
from textRepr import *
import time as wallTime
import sys

sys.path.append("C:\\Users\\Abhishek\\Desktop\\repositories\\libraryFunctions\\")

from parseAbqInpFileForNodalCoords import parseAbqFileForNodalCoordinates
from meltpoolCalculationsFromFE import findElementsContainingInterfaceSolidLiquid
from meltpoolCalculationsFromFE import findCubesEnclosingSolidLiquidInterfaceInAnElement
from meltpoolCalculationsFromFE import getCentroidValuesInterfaceCube

if len(sys.argv) != 6:
   text = """ Five arguments must be supplied : \
   \n\t (1) odbName with extension \
   \n\t (2) set name to extract the melt pool dimensions from \
   \n\t (3) variable Name (no quotations - NT11 for temperature) \
   \n\t (4) frameNoToExtractTheData \
   \n\t (5) liquidus temperature (should have same units as used in abq simulation)""" 
   raise RuntimeError(text)
    
start = wallTime.clock()    
# Settings
odbName = str(sys.argv[1])
setNameForMeltPoolExtraction = str(sys.argv[2])
variableName = str(sys.argv[3])
frameNo = int(sys.argv[4])
temperatureLiquidus = float(sys.argv[5])

# extract these from the odb
odb = openOdb(path = odbName)
stepNames = odb.steps.keys(0)
instanceName = odb.rootAssembly.instances.keys(0)
assembly = odb.rootAssembly
elType = assembly.instances[instanceName[0]].elements[0].type
if elType != 'DC3D8':
    raise RuntimeError('The mesh contains elements that are not D3CD8')

# Get different elsets
#elset_entireDomain = assembly.instances[instanceName[0]].elementSets['ENTIREGEOMETRY']
elset_meltPoolDimensions = assembly.instances[instanceName[0]].elementSets[setNameForMeltPoolExtraction.upper()]

# Extract data from the required frame
stepObject = odb.steps[stepNames[0]]
frameData = stepObject.frames[frameNo]
time = frameData.frameValue

print 'time: ', time

# Extract node IDs and nodal coordinates
nodalCoordObject = frameData.fieldOutputs['COORD']
nodalCoordFieldValues = nodalCoordObject.values
listCoordinatesNodal = [[],[],[]]
listNodeIDs = [int(iNodeField.nodeLabel) for iNodeField in nodalCoordFieldValues]
for iDimension in range(0,3):
    listCoordinatesNodal[iDimension] = [iNodeField.data[iDimension] for iNodeField in nodalCoordFieldValues]

# Extract nodal temperatures    
fieldVariableObject = frameData.fieldOutputs[variableName]
fieldVariableFieldValues = fieldVariableObject.values
listNodalTemperatures = [iFieldVariableField.data for iFieldVariableField in fieldVariableFieldValues]

# Extract element connectivity and required element sets
#listElementIDs_elsetEntireDomain = [iElementField.label for iElementField in elset_entireDomain.elements]
#listElementConnectivity_elsetEntireDomain = [iElementField.connectivity for iElementField in elset_entireDomain.elements]

listElementIDs = [iElementField.label for iElementField in elset_meltPoolDimensions.elements]
listElementConnectivity = [iElementField.connectivity for iElementField in elset_meltPoolDimensions.elements]

# Find elements that contain solid liquid interface
listElementsWithInterface, listInterfaceElementIndicesInElementSet = findElementsContainingInterfaceSolidLiquid(listElementIDs,
                                                                                                                listElementConnectivity,
                                                                                                                listNodalTemperatures,
                                                                                                                temperatureLiquidus)

if listElementsWithInterface == []:
    raise RuntimeError('There are zero elements on solid liquid interface')

# -------------------------------------- evaluate the centroid of each cube, temperature at the centroid and gradient at the centroid
listCoordinatesNaturalAtCentroidAll, listTemperaturesAtCentroidAll = [[],[],[]], []

for iElement in range(0,len(listElementsWithInterface)):

    idElement = listElementsWithInterface[iElement]
    
    indexElementInElSet = listInterfaceElementIndicesInElementSet[iElement]
    
    listNodesConnectedElement = listElementConnectivity[indexElementInElSet]
    
    listTemperaturesNodalElement = [listNodalTemperatures[idNode-1] for idNode in listNodesConnectedElement]
    
    listCoordinatesNodalElement = [[],[],[]]
    
    for iDimension in range(0,3):
    
        listCoordinatesNodalElement[iDimension] = [listCoordinatesNodal[iDimension][idNode-1] for idNode in listNodesConnectedElement]
    
    listNaturalCoodsInterfaceCubeNodes, temperatureNodesNewWithInterface = findCubesEnclosingSolidLiquidInterfaceInAnElement(listTemperaturesNodalElement,
                                                                                                                             temperatureLiquidus, 
                                                                                                                             maxDivisionsCube = 8)
                                                                                                                             
    for iCube in range(0, len(listNaturalCoodsInterfaceCubeNodes)):
        
        listCoordinatesCartesianAtCentroid, listCoordinatesNaturalAtCentroid, temperatureAtCentroid = getCentroidValuesInterfaceCube(listTemperaturesNodalElement,
                                                                                                                                     listCoordinatesNodalElement,
                                                                                                                                     listNaturalCoodsInterfaceCubeNodes[iCube])
    
        for iDimension in range(0,3):
        
            listCoordinatesNaturalAtCentroidAll[iDimension].append(listCoordinatesCartesianAtCentroid[iDimension])
            
        listTemperaturesAtCentroidAll.append(temperatureAtCentroid)
    
l_x = max(listCoordinatesNaturalAtCentroidAll[0])-min(listCoordinatesNaturalAtCentroidAll[0])
l_y = max(listCoordinatesNaturalAtCentroidAll[1])-min(listCoordinatesNaturalAtCentroidAll[1])
l_z = max(listCoordinatesNaturalAtCentroidAll[2])-min(listCoordinatesNaturalAtCentroidAll[2])

print('l_x,l_y,l_z : ', l_x, l_y, l_z)

# ---------------------------- write meltpool dimensions to a file
out_path = 'meltPoolDimensions.out'
with open(out_path, 'w') as file_out:    
    file_out.write("{0:25.10f}{1:25.10f}{2:25.10f}\n".format(l_x,l_y,l_z))        
file_out.close()

odb.close()
end = wallTime.clock()
print "Time Taken for extracting meltpool dimensions: ",(end-start), "seconds\n"