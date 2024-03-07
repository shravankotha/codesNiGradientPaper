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
   \n\t (2) abaqusInputFileName \
   \n\t (3) variable Name (no quotations - NT11 for temperature) \
   \n\t (4) frameNoToExtractTheData \
   \n\t (5) liquidus temperature (should have same units as used in abq simulation)""" 
   raise RuntimeError(text)
    
# Settings
odbName = str(sys.argv[1])
abaqusInputFileName = str(sys.argv[2])
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
elset_entireDomain = assembly.instances[instanceName[0]].elementSets['ENTIREGEOMETRY']
elset_deposit = assembly.instances[instanceName[0]].elementSets['SETDEPOSIT']
elset_basePlate = assembly.instances[instanceName[0]].elementSets['SETBASEPLATE']


# Extract data from the required frame
stepObject = odb.steps[stepNames[0]]
frameData = stepObject.frames[frameNo]
time = frameData.frameValue

# Extract node IDs and nodal coordinates
nodalCoordObject = frameData.fieldOutputs['COORD']
nodalCoordFieldValues = nodalCoordObject.values
listNodalCoordinates = [[],[],[]]
listNodeIDs = [int(iNodeField.nodeLabel) for iNodeField in nodalCcoordFieldValues]
for iDimension in range(0,3):
    listNodalCoordinates[iDimension] = [iNodeField.data[iDimension] for iNodeField in nodalCcoordFieldValues]
#listNodeIDs, listNodalCoordinates = parseAbqFileForNodalCoordinates(abaqusInputFileName)

# Extract nodal temperatures    
fieldVariableObject = frameData.fieldOutputs[variableName]
fieldVariableFieldValues = fieldVariableObject.values
listNodalTemperatures = [ifieldVariableField.data for iFieldVariableField in fieldVariableFieldValues]\

# Extract element connectivity and required element sets
listElementIDs_elsetEntireDomain = [iElementField.label for iElementField in elset_entireDomain.elements]
listElementConnectivity = [iElementField.connectivity for iElementField in elset_entireDomain.elements

# Find elements that contain solid liquid interface    
listElementsWithInterface = findElementsContainingInterfaceSolidLiquid(listElementIDs_elsetEntireDomain,
                                                                       listElementConnectivity,
                                                                       listNodalTemperatures,
                                                                       temperatureLiquidus)
                                                                       
# Find the min and max dimensions of the meltpool along X, Y and Z directions
minCoordAlongX, minCoordAlongY, minCoordAlongZ = 1E10,1E10,1E10
maxCoordAlongX, maxCoordAlongY, maxCoordAlongZ = -1E10,-1E10,-1E10

if listElementsWithInterface == []:
    raise RuntimeError('There zero elements on solid liquid interface')
    
for idElement in listElementsWithInterface:
    for idNode in listElementConnectivity[idElement-1]:
        if listCoordinatesNodal[idNode-1][0] < minCoordAlongX: 
            minCoordAlongX = listCoordinatesNodal[idNode-1][0]
            idElementWithMinX = idElement
        if listCoordinatesNodal[idNode-1][1] < minCoordAlongY:
            minCoordAlongX = listCoordinatesNodal[idNode-1][0]
            idElementWithMinY = idElement
        if listCoordinatesNodal[idNode-1][2] < minCoordAlongZ:
            minCoordAlongX = listCoordinatesNodal[idNode-1][0]
            idElementWithMinZ = idElement

        if listCoordinatesNodal[idNode-1][0] > maxCoordAlongX: 
            maxCoordAlongX = listCoordinatesNodal[idNode-1][0]
            idElementWithMaxX = idElement
        if listCoordinatesNodal[idNode-1][1] > maxCoordAlongY:
            maxCoordAlongX = listCoordinatesNodal[idNode-1][0]
            idElementWithMaxY = idElement
        if listCoordinatesNodal[idNode-1][2] > maxCoordAlongZ:
            maxCoordAlongX = listCoordinatesNodal[idNode-1][0]
            idElementWithMaxZ = idElement

# Find the maximum distance along X-axis



# -------------------------------------- evaluate the centroid of each cube, temperature at the centroid and gradient at the centroid

for idElement in [idElementWithMinX,idElementWithMaxX]:
    listNodesConnectedElement = listElementConnectivity[idElement-1]    
    listTemperaturesElement = [listNodalTemperatures[idNode-1] for idNode in listNodesConnectedElement]    
    listCoordinatesNodalElement = [[],[],[]]
    
    for iDimension in range(0,3):    
        listCoordinatesNodalElement[iDimension] = [listNodalCoordinates[iDimension][idNode-1] for idNode in listNodesConnectedElement]
    
    listNaturalCoodsInterfaceCubeNodes, temperatureNodesNewWithInterface = findCubesEnclosingSolidLiquidInterfaceInAnElement(listNodesConnectedElement,
                                                                                                                             listTemperaturesElement,
                                                                                                                             temperatureLiquidus)
    listCoordinatesCartesianAtCentroid, listCoordinatesNaturalAtCentroid, temperatureAtCentroid = getCentroidValuesInterfaceCube(listTemperaturesNodalElement,
                                                                                                                                 listCoordinatesNodalElement,
                                                                                                                                 listNaturalCoodsInterfaceCubeNodes)
    
    print('Coords, temp : ',listCoordinatesCartesianAtCentroid, temperatureAtCentroid)






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
end = wallTime.clock()
print "Time Taken for writing: ",(end-start), "seconds\n"