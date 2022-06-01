
import numpy as np

#objective function that we want to minimize
def objective_f(x, y):
    return (x**2 + y**2) / 400 - np.cos(x) * np.cos(y/(2**0.5)) + 1

#hyper parameters
crossover_rate = 0.99
mutation_rate = 0.3
maxgen = 1000
eps = 0.0
#
ax = -8
bx = 8
ay = -8
by = 8
M_population = 50

#initialize
def random_initialize():
    x_population = (bx - ax) * np.random.rand(M_population) + ax
    y_population = (by - ay) * np.random.rand(M_population) + ay
    return x_population, y_population


# tournament selection
def tournament_selection(obj_values, n=3):
    selected_ = np.random.randint(len(obj_values))
    for i in np.random.randint(0, len(obj_values), n-1):
        if obj_values[i] < obj_values[selected_]:
            selected_ = i
    return selected_


#arithmetic crossover
def arithmetic_crossover(P1, P2):
    C1 = np.copy(P1)
    C2 = np.copy(P2)
    if np.random.rand() < crossover_rate:
        alpha = np.random.rand()
        beta = 1 - alpha
        C1 = alpha * P1 + (1 - alpha) * P2
        C2 = beta * P1 + (1 - beta) * P2
    return C1, C2


def uniform_mutation(P):
    mutated = np.copy(P)
    if np.random.rand() < mutation_rate:
        mutated[0] = (bx - ax) * np.random.rand() + ax
        mutated[1] = (by - ay) * np.random.rand() + ay
    return mutated


def ga():
    x_best = ax
    y_best = ay
    obj_best = objective_f(ax, ay)
    x_population, y_population = random_initialize()
    generation = np.c_[x_population, y_population]
    obj_values = objective_f(generation[:, 0], generation[:, 1])
    sorted_idx = np.argsort(obj_values)
    for n_generation in range(maxgen):
        #save the elite one from the previous geneartion
        elite = generation[sorted_idx[0], :]
        elite_obj = obj_values[sorted_idx[0]]
        #check for improvement
        if elite_obj < obj_best:
            obj_best = elite_obj
            x_best = elite[0]
            y_best = elite[1]
        
        #do selection    
        selected = [tournament_selection(obj_values)
                    for i in range(M_population)]
        #do crossover
        crossovered_gen = np.zeros((M_population, 2))
        for i in range(0, M_population, 2):
            c1, c2 = arithmetic_crossover(generation[selected[i], :], generation[selected[i + 1], :])
            crossovered_gen[i, :] = c1
            crossovered_gen[i - 1, :] = c2
        
        #do mutation
        mutated_gen = np.zeros((M_population, 2))
        for i in range(M_population):
            mutated_gen[i, :] = uniform_mutation(crossovered_gen[i, :])
        #replacement
        generation = mutated_gen
        obj_values = objective_f(mutated_gen[:, 0], mutated_gen[:, 1])
        sorted_idx = np.argsort(obj_values)
        obj_values[sorted_idx[-1]] = elite_obj
        generation[sorted_idx[-1], :] = elite
        #log
        print("*"*20)
        print(f"generation number {n_generation + 1}, \nx_best = {x_best}, y_best = {x_best}, f(x_best, y_best) = {objective_f(x_best, y_best)}")
        #check for early convergence
        if obj_best <= eps:
            break
        
    return x_best, y_best

x_best, y_best = ga()
print("*"*20)
print(f"genetic algorithm has converged, final results:")
print(f'x_best = {x_best}, y_best = {x_best}, f(x_best, y_best) = {objective_f(x_best, y_best)}')
