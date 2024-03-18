### -------------------------------- Implement the problem 
import numpy as np
from pymoo.core.problem import ElementwiseProblem
import random
import time
import os
import shutil
from replaceAstringInAFile import replaceAstringInAFile

class MyProblem(ElementwiseProblem):

    def __init__(self):
        super().__init__(n_var = 2, n_obj = 2, n_ieq_constr = 0,
                         xl = np.array([0.1E-3,0.1E-3]),
                         xu = np.array([1.5E-3,1.5E-3]))

    def _evaluate(self, x, out, *args, **kwargs):
        pwd = os.getcwd()
        # create a new folder and copy reference simulation files
        source_dir_name = str(pwd) + '\\' + "RefFiles" + "\\"
        ref_inp_file_name = "beadAndHeatSourceInputs.inp"
        source_dir_path = str(source_dir_name)
        # change heat sources parameters in abaqus inp files
        random_number = int(time.time()*1000) # time since epoch in milliseconds
        dest_dir_name = 'Simulation_' + str(random_number)
        dest_dir_path = str(pwd) + '\\' + str(dest_dir_name) + "\\"
        if not os.path.exists(dest_dir_path):
            os.makedirs(dest_dir_path)
        
        shutil.copytree(source_dir_path, dest_dir_path, dirs_exist_ok=True)        
        nameFile = dest_dir_path + '\\' + str(ref_inp_file_name)
        
        beadWidth = 1.5E-3
        beadHeight = 0.6E-3
        laserSpotRadius = x[0]
        penetrationDepth = x[1]
        
        success = replaceAstringInAFile(nameFile,"__beadWidth__",str(beadWidth))
        success = replaceAstringInAFile(nameFile,"__beadHeight__",str(beadHeight))
        success = replaceAstringInAFile(nameFile,"__laserSpotRadius__",str(laserSpotRadius))
        success = replaceAstringInAFile(nameFile,"__penetrationDepth__",str(penetrationDepth))
        
        # run abaqus simulation
        command_to_execute = "abaqus job=singleTrack_Model.inp cpus=4 -inter"
        os.chdir(dest_dir_path)
        os.system(command_to_execute)
        
        # post-process to obtain melt-pool depth and width
        command_to_execute = "abaqus python C:\\Users\\Abhishek\\Desktop\\repositories\\codesNiGradientPaper\\getMeltpoolDimensionsFromAbaqusOutputs.py singleTrack_Model.odb setBasePlate NT11 4 1360"
        os.system(command_to_execute)
        
        [l_x,l_y,l_z] = np.loadtxt('meltPoolDimensions.out')        
        lengthSimulation = l_y
        widthSimulation = l_x
        depthSimulation = l_z
        print('laserSpotRadius, penetrationDepth, meltPoolWidth, meltPoolDepth : ', laserSpotRadius, penetrationDepth, widthSimulation, depthSimulation)
        os.chdir(pwd)
        
        # Evaluate the objective function
        # IN625
        widthTarget = 1.34E-3
        depthTarget = 0.29E-3
        # IN738
        #widthTarget = 1.49E-3
        #depthTarget = 0.33E-3
        
        f1 = (widthSimulation - widthTarget)**2
        f2 = (depthSimulation - depthTarget)**2
        
        out["F"] = [f1, f2]
        print ('-----------------------------------------------------------------------------------------------------------------------')
        
problem = MyProblem()

### -----------------------   Initialize an algorithm
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling

algorithm = NSGA2(pop_size = 4, n_offsprings = 2, sampling = FloatRandomSampling(),
    crossover = SBX(prob = 0.9, eta = 15), mutation = PM(eta = 20), eliminate_duplicates = True)
    
### ---------------------- Determine a termination criteria
from pymoo.termination import get_termination
termination = get_termination("n_gen", 2)

### ---------------------- Optimize
from pymoo.optimize import minimize
res = minimize(problem, algorithm, termination, seed = 1, save_history = True, verbose = True)    # problem, algorithm and termination are defined above

X = res.X
F = res.F
print('Optimal solution : ', X)
print('Objective function value : ', F)

### ---------------------- write the pareto front to file
out_path = 'paretoFront.out'
with open(out_path, 'w') as file_out:
    file_out.write("laserSpotRadius      penetrationDepth      objFunc_meltPoolWidth        objFunc_meltPoolDepth\n")
    for iRow in range(0, len(X[:,0])):
        file_out.write("{0:35.20f}{1:35.20f}{2:35.20f}{3:35.20f}\n".format(X[iRow,0], X[iRow,1], F[iRow,0], F[iRow,1]))       
file_out.close()

### --------------------- Visualize
import matplotlib.pyplot as plt
xl, xu = problem.bounds()
plt.figure(figsize=(7, 5))
plt.scatter(X[:, 0], X[:, 1], s=30, facecolors='none', edgecolors='r')
plt.xlim(xl[0], xu[0])
plt.ylim(xl[1], xu[1])
plt.title("Design Space")
plt.show()

plt.figure(figsize=(7, 5))
plt.scatter(F[:, 0], F[:, 1], s=30, facecolors='none', edgecolors='blue')
plt.title("Objective Space")
plt.show()