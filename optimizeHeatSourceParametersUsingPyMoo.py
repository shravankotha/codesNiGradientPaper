### -------------------------------- Implement the problem 
import numpy as np
from pymoo.core.problem import ElementwiseProblem
import random
import time
import os

class MyProblem(ElementwiseProblem):

    def __init__(self):
        super().__init__(n_var=2,
                         n_obj=2,
                         n_ieq_constr=0,
                         xl=np.array([-2,-2]),
                         xu=np.array([2,2]))

    def _evaluate(self, x, out, *args, **kwargs):
        pwd = os.getcwd()
        # create a new folder and copy reference simulation files
        source_dir_name = str(pwd) + '\\' + "RefFiles" + "\\"
        ref_inp_file_name = "beadAndHeatSourceInputs.inp"
        source_dir_path = str(source_dir_name)
        # change heat sources parameters in abaqus inp files        
        current_time = time.time()
        random.seed(current_time)
        random_number = random.randint(1, 1000000)        
        dest_dir_name = 'Simulation_' + str(random_number)
        dest_dir_path = str(pwd) + '\\' + str(dest_dir_name) + "\\"
        if not os.path.exists(dest_dir_path):
            os.makedirs(dest_dir_path)
        
        shutil.copytree(source_dir_path, dest_dir_path, dirs_exist_ok=True)        
        nameFile = dest_dir_path + '\\' + str(ref_inp_file_name)
        
        beadWidth = x[0]
        beadHeight = x[1]
        laserSpotRadius = x[2]
        penetrationDepth = x[3]
        
        success = replaceAstringInAFile(nameFile,"__beadWidth__",str(beadWidth))
        success = replaceAstringInAFile(nameFile,"__beadHeight__",str(beadHeight))
        success = replaceAstringInAFile(nameFile,"__laserSpotRadius__ ",str(laserSpotRadius))
        success = replaceAstringInAFile(nameFile,"__penetrationDepth__",str(penetrationDepth))
        
        # run abaqus simulation
        command_to_execute = "abaqus job=singleTrack_Model.inp cpus=4 -inter"
        os.chdir(dest_dir_path)
        os.system(command_to_execute)
        
        # post-process to obtain melt-pool depth and width
        
        
        # Evaluate the objective function        
        f1 = (depthSimulation1 - depthTarget1)**2
        f2 = (widthSimulation - widthTarget)**2
        out["F"] = [f1, f2]
        
problem = MyProblem()

### -----------------------   Initialize an algorithm
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling

algorithm = NSGA2(
    pop_size=40,
    n_offsprings=10,
    sampling=FloatRandomSampling(),
    crossover=SBX(prob=0.9, eta=15),
    mutation=PM(eta=20),
    eliminate_duplicates=True
)

### ---------------------- Determine a termination criteria
from pymoo.termination import get_termination
termination = get_termination("n_gen", 40)

### ---------------------- Optimize
from pymoo.optimize import minimize
res = minimize(problem,
               algorithm,
               termination,
               seed=1,
               save_history=True,
               verbose=True)

X = res.X
F = res.F
print('Optimal solution : ', X)
print('Objective function value : ', F)

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