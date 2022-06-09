
import numpy as np

C = np.array([[0, 8, 20, 35, 3, 15],
              [8, 0, 3, 5, 2, 1],
              [20, 3, 0, 20, 28, 9],
              [35, 5, 20, 0, 12, 22],
              [3, 2, 28, 12, 0, 10],
              [15, 1, 9, 22, 10, 0]], "d")

M_population = 15
num_cities = 6
crossover_rate = 0.99
mutation_rate = 0.1
maxgen = 100


class Edge:
    def __init__(self, c1, c2):
        self.c1 = c1
        self.c2 = c2
        self.C = C

    def __eq__(self, __o: int) -> bool:
        return self.c1 == __o or self.c2 == __o

    def __gt__(self, __o: object) -> bool:
        return self.C[self.c1 - 1, self.c2 - 1] > __o.C[__o.c1 - 1, __o.c2 - 1]

    def __lt__(self, __o: object) -> bool:
        return self.C[self.c1 - 1, self.c2 - 1] < __o.C[__o.c1 - 1, __o.c2 - 1]

    def __ge__(self, __o: object) -> bool:
        return self.C[self.c1 - 1, self.c2 - 1] >= __o.C[__o.c1 - 1, __o.c2 - 1]

    def __le__(self, __o: object) -> bool:
        return self.C[self.c1 - 1, self.c2 - 1] <= __o.C[__o.c1 - 1, __o.c2 - 1]

    def __add__(self, __o: object):
        if self.c1 == __o.c1:
            if self.c2 == __o.c2:
                return Edge(self.c1, self.c2)
            else:
                return Edge(self.c1, 0)
        elif self.c1 == __o.c2:
            if self.c2 == __o.c1:
                return Edge(self.c1, self.c2)
            else:
                return Edge(self.c1, 0)
        elif self.c2 == __o.c2:
            return Edge(0, self.c2)
        elif self.c2 == __o.c1:
            return Edge(0, self.c2)
        else:
            return Edge(0, 0)

    def __str__(self):
        return f'{self.c1}->{self.c2}'


class EdgePath:
    def __init__(self, path):
        self.path = path
        self.update()

    def update(self):
        self.edges = []
        if not len(self.path) == 0:
            for i in range(len(self.path) - 1):
                self.edges.append(Edge(self.path[i], self.path[i + 1]))
            self.edges.append(Edge(self.path[-1], self.path[0]))

    def __add__(self, __o: object):
        com = []
        for e1 in self.edges:
            for e2 in __o.edges:
                if e1 == e2:
                    if not e1 in com:
                        com.append(e1)
        return com

    def to_array(self):
        return self.path

    def __str__(self) -> str:
        s = "|"
        for e in self.edges:
            s += f'{e})('
        s += '|'
        return s


class Population:

    def __init__(self, M_population=15, n_genes=6):

        self.C = C
        self.M_population = M_population
        self.n_genes = n_genes

        self.population = self.init_path_code()
        self.update_obj_values()
        self.update_edges()

    def init_path_code(self):
        population = []
        for _ in range(self.M_population):
            population.append(np.random.permutation(self.n_genes) + 1)
        return np.array(population, 'i')

    def update_edges(self):
        self.edge_path_list = []
        for x in self.population:
            self.edge_path_list.append(EdgePath(x))

    def objective_func(self, x):
        z = 0.0
        for i in range(len(x) - 1):
            z += self.C[int(x[i] - 1), int(x[i + 1] - 1)]
        z += self.C[int(x[-1] - 1), int(x[0] - 1)]
        return z

    def update_obj_values(self):
        self.obj_values = [self.objective_func(x) for x in self.population]

    def update(self):
        self.update_obj_values()
        self.update_edges()


def find_edges(edge_list, num):
    b = []
    for e in edge_list:
        if e == num:
            b.append(e)
    if b:
        bs = np.sort(b)
        return bs[0].c1 if bs[0].c2 == num else bs[0].c2
    else:
        return None


def edge_recombination_crossover(x: EdgePath, y: EdgePath) -> EdgePath:
    if np.random.rand() < crossover_rate:
        p = [0] * num_cities
        edges = []
        for i in range(1, num_cities + 1):
            for j in range(num_cities):
                if x.edges[j] == i:
                    if not x.edges[j] in edges:
                        p[i - 1] += 1
                        edges.append(x.edges[j])

                if y.edges[j] == i:
                    if not y.edges[j] in edges:
                        p[i - 1] += 1
                        edges.append(y.edges[j])

        a_sort = np.argsort(p) + 1
        c = x + y
        gene = []
        current = a_sort[0]
        while True:
            if not current in gene:
                gene.append(current)

            if(len(gene) == len(x.path)):
                break

            b = find_edges(c, current)
            if b and not b in gene:
                current = b
            else:
                b = find_edges(edges, current)
                if b and not b in gene:
                    current = b
                else:
                    b = current
                    while b in gene:
                        b = np.random.randint(1, 7)
                    current = b
        return gene
    else:
        return x.path


def tournament_selection(obj_values, n=3):
    selected_ = np.random.randint(len(obj_values))
    for i in np.random.randint(0, len(obj_values), n-1):
        if obj_values[i] < obj_values[selected_]:
            selected_ = i
    return selected_


def simple_inversion_mutation(p):
    if np.random.rand() < mutation_rate:
        new_p = np.copy(p)
        r1 = np.random.randint(0, len(p) - 2)
        r2 = np.random.randint(2, len(p))
        c2 = max(r1, r2)
        c1 = min(r1, r2)
        if c2 == c1:
            c2 = c1 + 2
        x = new_p[c1:c2]
        new_p[c1:c2] = x[::-1]
        return new_p
    else:
        return p


def ga():
    generation = Population(M_population=M_population, n_genes=num_cities)
    obj_values = generation.obj_values
    sorted_idx = np.argsort(obj_values)
    for n_generation in range(maxgen):
        # save the elite one from the previous geneartion
        elite = generation.population[sorted_idx[0], :]
        elite_obj = obj_values[sorted_idx[0]]

        # do selection
        selected = [tournament_selection(obj_values)
                    for i in range(M_population)]
        # do crossover
        crossovered_gen = np.zeros((M_population, num_cities), "i")
        for i in range(M_population - 1):

            ep = edge_recombination_crossover(
                generation.edge_path_list[selected[i]], generation.edge_path_list[selected[i + 1]])
            crossovered_gen[i, :] = ep
        crossovered_gen[-1, :] = elite
        # do mutation
        mutated_gen = np.zeros((M_population, 6), 'i')
        for i in range(M_population):
            mutated_gen[i, :] = simple_inversion_mutation(
                crossovered_gen[i, :])
        # replacement
        generation.population = np.array(mutated_gen, "i")
        generation.update()
        obj_values = generation.obj_values
        sorted_idx = np.argsort(obj_values)
        generation.population[sorted_idx[-1], :] = elite

        generation.obj_values[sorted_idx[-1]] = elite_obj
        # log

    return elite, generation


elite, generation = ga()
print("*"*20)
print(f"genetic algorithm has converged, final results:")
print(elite)

print(EdgePath(elite))
