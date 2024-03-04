### -------------------------------- Implement the problem 
import numpy as np
from pymoo.core.problem import ElementwiseProblem

class MyProblem(ElementwiseProblem):

    def __init__(self):
        super().__init__(n_var=2,
                         n_obj=2,
                         n_ieq_constr=0,
                         xl=np.array([-2,-2]),
                         xu=np.array([2,2]))

    def _evaluate(self, x, out, *args, **kwargs):
        # create a new folder and copy reference simulation files
        
        # change heat sources parameters in abaqus inp files
        
        
        # run abaqus simulation
        
        
        # post-process to obtain melt-pool depth and width
        
        
        # Evaluate the objective function
        
        f1 = (depthSimulation - depthTarget)**2
        f2 = (widthSimulation - widthTarget)**2
        
        #g1 = 2*(x[0]-0.1) * (x[0]-0.9) / 0.18
        #g2 = - 20*(x[0]-0.4) * (x[0]-0.6) / 4.8

        out["F"] = [f1, f2]
        #out["G"] = [g1, g2]
        
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