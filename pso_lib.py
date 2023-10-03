import random
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy

random.seed(1905)

class optimiser:
    def __init__(self, dimensions)-> None:
        self.dimensions = dimensions
    def fitness_funtion(self, position):
        pass

class ackley(optimiser):
    def fitness_function(self, position):
        first_summation = 0
        sec_summation = 0
        for i in range(self.dimensions):
            first_summation += position[i]**2
            sec_summation += np.cos(2 * np.pi * position[i])
        first_term = -20 * np.exp(-0.2 * math.sqrt((1/self.dimensions) * first_summation))
        sec_term = -np.exp((1/self.dimensions) * sec_summation) + 20 + np.exp(1)
        result = first_term + sec_term
        return result
    class rosenbrock(optimiser):
        def fitness_function(self, position):
            n = self.dimensions
            sum_val = 0.0
            for i in range(n - 1):
                sum_val += 100 * (position[i + 1] - position[i]**2)**2 + (1 - position[i])**2
            return sum_val

class PSO:
    def __init__(self, size, dimensions, topology, adjustment) -> None:
        self.opt = ackley(dimensions)
        self.size = size
        self.dimensions = dimensions
        self.topology = topology
        self.adjustment = adjustment
        self.swarm = swarm(-30, 30, size, dimensions, topology, adjustment, self.opt)
        self.gbest_list = []
        self.distance_mass_centre = []
        self.standard_devs = []
        self.velocity_vec_length = []

    def standard_PSO(self, iterations):
        for i in range(iterations):
            for ind in self.swarm.pop:
                ind.update_pos(i, iterations, 0)
            self.track_stats()
                
    def sync_PSO(self, iterations):
        for i in range(iterations):
            for ind in self.swarm.pop:
                ind.update_pbest([ind.position, ind.fit])
                ind.update_gbest()
            for ind in self.swarm.pop:
                ind.update_pos(i, iterations, 1)
            self.track_stats()
            if i % 10 == 0:
                print(self.swarm.pop[0].pbest[1])
        self.reset()

    def printer(self):
        for n in self.swarm.pop:
            print(n.position, n.neighbours)

    def reset(self):
        self.swarm = swarm(-30, 30, self.size, self.dimensions, self.topology, self.adjustment, self.opt)
    
    def clear_lists(self):
        self.gbest_list.clear()
        self.velocity_vec_length.clear()
        self.standard_devs.clear()
        self.distance_mass_centre.clear()
    def track_stats(self):
        self.gbest_list.append(self.swarm.find_true_global_best()) #works
        self.distance_mass_centre.append(self.swarm.get_dist()) # works
        self.standard_devs.append(self.swarm.get_standard_devs())
        self.velocity_vec_length.append(self.swarm.get_mean_len())


class swarm:
    def __init__(self, l_boundary, u_boundary, size, dimensions, topology, adjustment, opt):
        self.size = size
        self.u_boundary = u_boundary
        self.l_boundary = l_boundary
        self.dimensions = dimensions
        self.topology = topology,
        self.adjustment = adjustment
        self.opt = opt
        self.pop = []
        self.make_swarm()

    def make_swarm(self):
        for _ in range(self.size):
            self.pop.append(particle(self.dimensions, self.l_boundary, self.u_boundary, self.adjustment,self.opt))
        potential_neighbours = copy.deepcopy(self.pop)
        self.assign_neighbours(potential_neighbours)

    def get_ind(self, index):
        return self.pop[index]
    
    def get_mean_len(self):
        total_len = 0
        for ind in self.pop:
            total_len += np.linalg.norm(ind.velocity)
        return total_len/self.size
    
    def find_true_global_best(self):
        best = self.pop[0]
        for ind in self.pop:
            if ind.gbest[1] < best.gbest[1]:
                best = ind
        return best.gbest[1]

    def get_rand_ind(self):
        return self.get_ind(random.randint(0, self.size - 1))
    
    def assign_neighbours(self, potential_neighbours):
        if self.topology[0] ==  "lbest":
            for i in range(1, self.size - 1):
                self.pop[i].add_neighbours([self.get_ind(i-1), self.get_ind(i+1)])
            self.pop[0].add_neighbours([self.get_ind(1), self.get_ind(-1)])
            self.pop[-1].add_neighbours([self.get_ind(0), self.get_ind(-2)])

        elif self.topology[0] ==  "gbest":
            for ind in self.pop:
                ind.add_neighbours(potential_neighbours)

        elif self.topology[0] ==  "star":
            centre = self.get_rand_ind()
            centre.add_neighbours(potential_neighbours)
            for ind in self.pop:
                if ind is not centre:
                    ind.add_neighbours([centre])

        elif self.topology[0] ==  "rand":
           num = random.randint(1, self.size - 1)
           for ind in self.pop:
               for _ in range(num):
                   n = self.get_rand_ind()
                   ind.add_neighbours([n])
                   n.add_neighbours([ind])
        else:
            raise Exception("Not a valid topology.")

    def find_mean_pos(self):
        mean_pos = [0] * self.dimensions
        for i in range(self.dimensions):
            for ind in self.pop:
                mean_pos[i] += ind.position[i]
            mean_pos[i] = mean_pos[i]/self.size
        return mean_pos
 
    def euclidean_distance(self, point1):
        point2 = [0]*self.dimensions
        squared_diff = 0.0
    
    # Calculate the sum of squared differences for each dimension
        for i in range(len(point1)):
            squared_diff += (point1[i] - point2[i]) ** 2
        
        # Take the square root of the sum of squared differences to get the distance
        distance = squared_diff ** 0.5
        #print(distance)
        return distance
    
    def get_dist(self):
        centre_of_mass = self.find_mean_pos()
        return self.euclidean_distance(centre_of_mass)

    def get_standard_devs(self):
        distances = []
        for ind in self.pop:
            distances.append(self.euclidean_distance(ind.position))
        mean_pos = self.find_mean_pos()
        mean_dist = self.euclidean_distance(mean_pos)

        std_deviations = []
        for i in range(self.size):
            std_deviations.append((distances[i] - mean_dist)**2)
        variance = 0
        for i in range(len(std_deviations)):
            variance += std_deviations[i]
        variance /= self.size
        std_deviation = math.sqrt(variance)
        return std_deviation

class particle:
    def __init__(self, dimensions, l_boundary, u_boundary,adjustment, optimiser):
        self.position = self.get_rand_num(l_boundary + 1, u_boundary - 1, dimensions)
        self.dimensions = dimensions
        self.velocity = self.get_rand_num(1, 1.02, dimensions)
        self.optimiser = optimiser
        self.adjustment = adjustment
        self.fit = self.optimiser.fitness_function(self.position)
        self.neighbours = []
        self.pbest = [self.position, self.fit]
        self.gbest = [self.position, self.fit]
        

    def calculate_fit(self):
        self.fit = self.optimiser.fitness_function(self.position)
    
    def get_pbest(self):
        return self.pbest
    def get_rand_num(self, lower, upper, number):
        rand_list = []
        for _ in range(number):
            rand_list.append(random.uniform(lower, upper))
        return rand_list
    
    def add_neighbours(self, new_neighbours):
        for prospect in new_neighbours:
            if prospect in self.neighbours or prospect == self:
                return
        self.neighbours += new_neighbours
    
    def update_pos(self, curr_it, total_it, sync):
        if self.adjustment == "constriction":
            self.constrict_update(sync)
        elif self.adjustment == "nonlinear inertia":
            self.nonlin_inertia_update(curr_it, total_it, sync)
        elif self.adjustment == "linear inertia":
            self.lin_inertia_update(curr_it, total_it, sync)
        else:
            raise Exception("Not a valid update method.")

    def update_pbest(self, new_data):
       if new_data[1] < self.pbest[1]:
           self.pbest = new_data
            

    def update_gbest(self):
        gbest = self.pbest
        for ind in self.neighbours:
            if ind.pbest[1] < gbest[1]:
                gbest = ind.pbest
        self.gbest = gbest
        
    # We do not need this one, but might be interesting anyway :-)
    def lin_inertia_update(self, curr_it, total_it, sync):
        c12 = 2.05
        e1 = random.random()
        e2 = random.random()
        pbest_c = self.get_pbest()[0]
        w_min = 0.4
        w_max = 0.9
        w = w_max - ((w_max - w_min) * curr_it/total_it)
        new_position = [10] * self.dimensions
        for d in range (self.dimensions):
            self.velocity[d] = w * self.velocity[d] + c12*e1*(pbest_c[d]- self.position[d]) + c12*e2*(self.gbest[0][d]- self.position[d])
            new_position[d] = self.position[d] + self.velocity[d]
        self.position = new_position
        self.calculate_fit()
        if sync == 0:
            self.update_pbest([new_position, self.fit] )
            self.update_gbest()
    # This is the variation using non-linear inertia
    def nonlin_inertia_update(self, curr_it, total_it, sync):
        c12 = 2.05
        e1 = random.random()
        e2 = random.random()
        pbest_c = self.get_pbest()[0]
        w = self.calc_omega(curr_it, total_it)
        new_position = [10] * self.dimensions
        for d in range (self.dimensions):
            self.velocity[d] = w * self.velocity[d] + c12*e1*(pbest_c[d]- self.position[d]) + c12*e2*(self.gbest[0][d]- self.position[d])
            new_position[d] = self.position[d] + self.velocity[d]
        self.position = new_position
        self.calculate_fit()
        if sync == 0:
            self.update_pbest([new_position, self.fit] )
            self.update_gbest()
        
    
    
    def calc_omega(self, t, T):
        n = 1
        w_min = 0.8 #good
        w_max = -0.1 #good
        w = ((pow((T - t)/T, n ))*(w_min - w_max)) + w_max
        return w
    
    # This is the one for standard PSO. All numbers are from the paper.
    def constrict_update(self, sync):
        xi = 0.72984
        c12 = 2.05
        e1 = random.random()
        e2 = random.random()
        pbest_c = self.get_pbest()[0]
        new_position = [10] * self.dimensions
        for d in range (self.dimensions):
            self.velocity[d] = xi*(self.velocity[d] + c12*e1*(pbest_c[d]- self.position[d]) + c12*e2*(self.gbest[0][d]- self.position[d]))
            new_position[d] = self.position[d] + self.velocity[d]
        self.position = new_position
        self.calculate_fit()
        if sync == 0:
            self.update_pbest([new_position, self.fit] )
            self.update_gbest()

def graph(data, reps, msg):
    data1 = []
    data2 = []
    data3 = []
    data4 = []
    colors = cm.get_cmap('tab10', reps)
    legend_labels = ["run 1", "run 2", "run 3", "run 4", "run 5", "run 6", "run 7", "run 8", "run 9", "run 10"]
    for i in range(len(data)):
        data1.append(data[i][0])
        data2.append(data[i][1])
        data3.append(data[i][2])
        data4.append(data[i][3])
    
    plt.figure(figsize=(20, 8))
    co_it = 0
    
    plt.subplot(2, 2, 1)
    for listy, label in zip(data1, legend_labels):
        plt.plot(range(1, 1001), listy, color=colors(co_it), label=label)
        co_it += 1
    plt.xlabel('Iterations')
    plt.ylabel('Fitness')
    plt.title('Global best over 1000 iterations')
    plt.grid(alpha=0.5, linestyle='--')
    plt.legend( fontsize="5")
    
    co_it = 0
    plt.subplot(2, 2, 2)
    for listy, label in zip(data2, legend_labels):
        plt.plot(range(1, 1001), listy, color=colors(co_it), label=label)
        co_it += 1
    plt.xlabel('Iterations')
    plt.ylabel('Distance')
    plt.title('Distance of swarm mean position to origin')
    plt.grid(alpha=0.5, linestyle='--')
    plt.legend(fontsize="5")
    co_it = 0
    plt.subplot(2, 2, 3)
    for listy, label in zip(data3, legend_labels):
        plt.plot(range(1, 1001), listy, color=colors(co_it), label=label)
        co_it += 1
    plt.xlabel('Iterations')
    plt.ylabel('Standard deviation')
    plt.title('Standard deviation of particles to mean position')
    plt.grid(alpha=0.5, linestyle='--')
    plt.legend(fontsize="5")
    co_it = 0
    plt.subplot(2, 2, 4)
    for listy, label in zip(data4, legend_labels):
        plt.plot(range(1, 1001), listy, color=colors(co_it), label=label)
        co_it += 1
    plt.ylabel('Mean vector length')
    plt.xlabel('Iterations')
    plt.title('Mean velocity of swarm')
    plt.grid(alpha=0.5, linestyle='--')
    plt.legend(fontsize="5")

    plt.subplots_adjust(hspace=0.5)
    plt.savefig('metrics_std_pso{}.png'.format(msg))
    plt.close()
