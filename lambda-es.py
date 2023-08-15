import numpy as np
from numpy import sin, cos, pi
from scipy.optimize import root
from es import CMAES

dims_0 = (29, 49, 39, 40, 500, 0)
npts = 50
epsilon_0 = 1
TR_cutoff_0 = 42
guess_0 = np.array([-0.37, 0.93, 0.47, -0.88]) #almost 400mm test stand 


def douglas_peucker(points, epsilon):
    if len(points) <= 2:
        return points
    
    dmax = 0
    index = 0
    
    for i in range(1, len(points) - 1):
        d = np.abs(np.cross(points[-1] - points[0], points[i] - points[0])) / np.linalg.norm(points[-1] - points[0])
        if d > dmax:
            index = i
            dmax = d
    
    if dmax > epsilon:
        recursive_result1 = douglas_peucker(points[:index + 1], epsilon)
        recursive_result2 = douglas_peucker(points[index:], epsilon)
        return np.concatenate((recursive_result1[:-1], recursive_result2))
    else:
        return np.array([points[0], points[-1]])
    
def get_endpt(dims, theta1, guess):
    L1, L2, L3, L4, L5, theta5 = dims
    sin5 = sin(theta5)
    cos5 = cos(theta5)

    sin1 = sin(theta1)
    cos1 = cos(theta1) 
    def func_sincos(vars):
        sin2, cos2, sin3, cos3 = vars
        return [
            L1*sin1 + L2*sin2 + L3*sin3,
            L1*cos1 + L2*cos2 + L3*cos3 - L4, 
            sin2**2 + cos2**2 - 1,
            sin3**2 + cos3**2 - 1,        
        ]
    def jac_sincos(vars):
        sin2, cos2, sin3, cos3 = vars
        return [
            [L2, 0, L3, 0],
            [0, L2, 0, L3],
            [2*sin2, 2*cos2, 0, 0],
            [0, 0, 2*sin3, 2*cos3]
        ]
    
    sol = root(func_sincos, jac=jac_sincos, x0=guess, method='hybr').x
    sin2, cos2, sin3, cos3 = sol
    guess = sol

    p5 = L1*np.array([cos1, sin1]) + L2*np.array([cos2, sin2]) + L5*np.array([cos2*cos5-sin2*sin5, sin2*cos5+cos2*sin5])
    return p5, sol

def cost(dims, guess, npts=50, TR_cutoff=TR_cutoff_0, epsilon=epsilon_0):
    theta1_range = np.linspace(0, 2*pi, npts)
    path = []

    for theta1 in theta1_range:
        p5, sol = get_endpt(dims, theta1, guess)
        guess = sol
        path.append(p5)

    start_i = int(npts/2)
    end_i = int(npts/2)

    last_p5 = path[start_i]
    while start_i > 0:
        p5 = path[start_i]
        if np.linalg.norm(p5 - last_p5) > TR_cutoff:
            break
        last_p5 = p5
        start_i -= 1

    last_p5 = path[end_i]
    while end_i < npts:
        p5 = path[end_i]
        if np.linalg.norm(p5 - last_p5) > TR_cutoff:
            break
        last_p5 = p5
        end_i += 1

    path = np.array(path[start_i:end_i])
    path_simp = douglas_peucker(path, epsilon)

    max_length = 0
    for i in range(len(path_simp) - 1):
        segment_length = np.linalg.norm(path_simp[i] - path_simp[i + 1])
        max_length = max(segment_length, max_length)
    
    return -max_length


print("initial cost:", cost(dims_0, guess_0))


NPARAMS = 3        # make this a 100-dimensinal problem.
NPOPULATION = 20    # use population size of 101.
MAX_ITERATION = 10 # run each solver for 5000 generations.

# defines a function to use solver to solve fit_func
def test_solver(solver):
  history = []
  for j in range(MAX_ITERATION):
    solutions = solver.ask()
    # if j == 0:
    #   solutions = np.full((NPOPULATION, NPARAMS), guess)

    # print(solutions)

    fitness_list = np.zeros(solver.popsize)
    for i in range(solver.popsize):
      fitness_list[i] = cost(solutions[i])

    solver.tell(fitness_list)
    result = solver.result() # first element is the best solution, second element is the best fitness
    history.append(result[1])
    if (j+1) % 1 == 0:
      print("fitness at iteration", (j+1), result[1], result[0], solver.rms_stdev())
  print("local optimum discovered by solver:\n", result[0])
  print("fitness score at this local optimum:", result[1])
  return history

# defines CMA-ES algorithm solver
cmaes = CMAES(NPARAMS,
  x0=guess_0,
  popsize=NPOPULATION,
  weight_decay=0.0,
  sigma_init = 0.01
)



cma_history = test_solver(cmaes)

print(cma_history)