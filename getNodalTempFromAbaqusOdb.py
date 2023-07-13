'''
This code extracts connectivity information and the nodal temperatures for all the abaqus steps at all time increments
'''
from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
from textRepr import *
import time as wallTime

if len(sys.argv) != 4:
   text = """Three command line arguments are expected: \
            \n\t(1) odbName with extension \
            \n\t(2) output directory path \
            \n\t(3) frame increment Number"""
   raise RuntimeError(text)
    
# Settings
odbName = str(sys.argv[1])
outputPath = sys.argv[2]
frameNoInc = int(sys.argv[3])             # useful if no of frames are too large - first and last frames are always exported

# extract these from the odb
odb = openOdb(path = odbName)
start = wallTime.clock()
for iStep in range(0,len(odb.steps)):
    time_start = wallTime.clock()
    print('-----------------------------------------------------------')    
    stepName = odb.steps.keys(0)
    stepObject = odb.steps[stepName[iStep]]
    print('Step : ' + str(iStep+1) + '/' + str(len(odb.steps)) + ' -- ' + str(stepName[iStep]))
    instanceName = odb.rootAssembly.instances.keys(0) # 0 indicates first instance
    assembly = odb.rootAssembly
    elType = assembly.instances[instanceName[0]].elements[0].type
  
    if elType != 'DC3D8':
        print 'Element Type: ',elType 
        raise RuntimeError('This code works only for DC3D8 elements')
    
    if iStep == 0:        
        file_elemConn = open(outputPath + '/elementConnectivity.dat','w')        
    file_temperatures = open(outputPath + '/time_temps_' + str(stepName[iStep]) + '.dat','w')

    totalNoFrames = len(stepObject.frames)           
    frameNosToExport = range(0, totalNoFrames-1, frameNoInc)
    frameNosToExport.append(totalNoFrames-1)  # this ensures that the last time step is exported
    noOfFramesToExport = len(frameNosToExport)    

    for frameNumber in range(1, len(frameNosToExport)):
        frameData = stepObject.frames[frameNosToExport[frameNumber]]
        time = frameData.frameValue
        print(' frame Number: ' + str(frameNumber) + ' time: ' + str(time))
        if iStep == 0 and frameNumber == 1:            
            for name, instance in assembly.instances.items():
                for element in instance.elements:
                    for nodeNum in element.connectivity:
                        file_elemConn.write('%7d'%nodeNum)
                    file_elemConn.write('\n')
            
        # Write temperatures at all times for this step 
        file_temperatures.write('time : ' + str(time) + '\n')    
        temperatureObject = frameData.fieldOutputs['NT11']        
        temperatureFieldValues = temperatureObject.values
        for iTempField in temperatureFieldValues:
            file_temperatures.write('%14.5F'%iTempField.data)
            file_temperatures.write('\n')
        
    file_temperatures.close()
    if iStep == 0:    
       file_elemConn.close()
    time_end = wallTime.clock()
    print('  -- time elapsed for the step : ' + str(time_end-time_start) + ' seconds')
odb.close()
end = wallTime.clock()
print "Total Time Taken: ",(end-start), "seconds\n"