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
from meltpoolCalculationsFromFE import mapIdElementToiElementAbaqus

if len(sys.argv) != 5:
   text = """ Four arguments must be supplied : \
   \n\t (1) odbName with extension \
   \n\t (2) variable Name (no quotations - NT11 for temperature) \
   \n\t (3) frameNoToExtractTheData \
   \n\t (4) liquidus temperature (should have same units as used in abq simulation)""" 
   raise RuntimeError(text)
    
# Settings
odbName = str(sys.argv[1])
variableName = str(sys.argv[2])
frameNo = int(sys.argv[3])
temperatureLiquidus = float(sys.argv[4])

# extract these from the odb
odb = openOdb(path = odbName)
stepNames = odb.steps.keys(0)
instanceName = odb.rootAssembly.instances.keys(0)
assembly = odb.rootAssembly
elType = assembly.instances[instanceName[0]].elements[0].type
if elType != 'DC3D8':
    raise RuntimeError('The mesh contains elements that are not D3CD8')

# Get different elsets
elset_entireDomain = assembly.instances[instanceName[0]].elementSets['ENTIREGEOMETRY']
elset_deposit = assembly.instances[instanceName[0]].elementSets['SETDEPOSIT']
elset_basePlate = assembly.instances[instanceName[0]].elementSets['SETBASEPLATE']


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
#listNodeIDs, listCoordinatesNodal = parseAbqFileForNodalCoordinates(abaqusInputFileName)

# Extract nodal temperatures    
fieldVariableObject = frameData.fieldOutputs[variableName]
fieldVariableFieldValues = fieldVariableObject.values
listNodalTemperatures = [iFieldVariableField.data for iFieldVariableField in fieldVariableFieldValues]

# Extract element connectivity and required element sets
listElementIDs_elsetEntireDomain = [iElementField.label for iElementField in elset_entireDomain.elements]
listElementConnectivity = [iElementField.connectivity for iElementField in elset_entireDomain.elements]

# Find elements that contain solid liquid interface
listElementsWithInterface = findElementsContainingInterfaceSolidLiquid(listElementIDs_elsetEntireDomain,
                                                                       listElementConnectivity,
                                                                       listNodalTemperatures,
                                                                       temperatureLiquidus)
                                                          
mapIdElementToiElement = mapIdElementToiElementAbaqus(listElementIDs_elsetEntireDomain)

# Find the min and max dimensions of the meltpool along X, Y and Z directions
minCoordAlongX, minCoordAlongY, minCoordAlongZ = 1E10,1E10,1E10
maxCoordAlongX, maxCoordAlongY, maxCoordAlongZ = -1E10,-1E10,-1E10

if listElementsWithInterface == []:
    raise RuntimeError('There are zero elements on solid liquid interface')

# -------------------------------------- evaluate the centroid of each cube, temperature at the centroid and gradient at the centroid
listCoordinatesNaturalAtCentroidAll, listTemperaturesAtCentroidAll = [[],[],[]], []

for idElement in listElementsWithInterface:

    iElement = mapIdElementToiElement[idElement-1]
    
    listNodesConnectedElement = listElementConnectivity[iElement]
    
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

out_path = 'meltPoolDimensions.out'

with open(out_path, 'w') as file_out:
    
    file_out.write("{0:25.10f}{1:25.10f}{2:25.10f}\n".format(l_x,l_y,l_z)) 
       
file_out.close()

#elementIDsObject = elset_deposit.elements[0]
#listElementIDs_elsetDeposit = []
#for iElement in range(len(elementIDsObject)):
#    listElementIDs_elsetDeposit.append(elementIDsObject[iElement].label)
#
#elementIDsObject = elset_basePlate.elements[0]
#listElementIDs_elsetDeposit = []
#for iElement in range(len(elementIDsObject)):
#    listElementIDs_elsetDeposit.append(elementIDsObject[iElement].label)




# Write melt-pool dimensions to file
#outputFileName = 'meltPoolDimensions.dat'
#file_output = open(outputFileName,'w')
#baseString = '  Time     Length     Width       Depth\n'
#file_output.write(baseString)
#file_output.write('%20.8E\t%20.8E\t%20.8E\t\n',data_to_write)   
#file_output.close()



odb.close()
#end = wallTime.clock()
#print "Time Taken for writing: ",(end-start), "seconds\n"