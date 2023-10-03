import pso_lib as PSO
import copy

# to do:
#make sure gbest does not include the individual!
def run_exp(size, filename, topology, adjustment, reps, gens):
    standardPSO = PSO.PSO(size, 10, topology, adjustment)
    data = []
    for _ in range(10):
        standardPSO.standard_PSO(gens)
        data.append(copy.deepcopy([standardPSO.gbest_list, standardPSO.distance_mass_centre, standardPSO.standard_devs, standardPSO.velocity_vec_length]))
        standardPSO.clear_lists()
        standardPSO.reset()
    PSO.graph(data, reps, filename)

i = 1
run_exp(50, "", "lbest", "constriction", 10, 1000)
print(i)
i += 1
run_exp(50, "_topo_gbest", "gbest", "constriction", 10, 1000)
print(i)
i += 1
run_exp(50, "_topo_rand", "rand", "constriction", 10, 1000)
print(i)
i += 1
run_exp(50, "_topo_star", "star", "constriction", 10, 1000)
print(i)
i += 1
run_exp(50, "_weight_adjust", "lbest", "nonlinear inertia", 10, 1000)
print(i)
i += 1
run_exp(20, "20", "lbest", "constriction", 10, 1000)
print(i)
i += 1
run_exp(100, "100", "lbest", "constriction", 10, 1000)
print(i)
i += 1
run_exp(200, "200", "lbest", "constriction", 10, 1000)
print(i)
i += 1
#for 10) analysis
run_exp(20, "_weight_adjust_20", "lbest", "nonlinear inertia", 10, 1000)
print(i)
i += 1

run_exp(20, "_weight_adjust_topo_gbest_20", "gbest", "nonlinear inertia", 10, 1000)
print(i)
i += 1
run_exp(20, "_topo_gbest_20", "gbest", "constriction", 10, 1000)
print(i)
i += 1

run_exp(20, "_weight_adjust_topo_star_20", "star", "nonlinear inertia", 10, 1000)
print(i)
i += 1
run_exp(20, "_topo_star_20", "star", "constriction", 10, 1000)
print(i)
i += 1

run_exp(20, "_weight_adjust_topo_random_20", "rand", "nonlinear inertia", 10, 1000)
print(i)
i += 1
run_exp(20, "_topo_random_20", "rand", "constriction", 10, 1000)
print(i)
i += 1


run_exp(100, "_weight_adjust_100", "lbest", "nonlinear inertia", 10, 1000)
print(i)
i += 1

run_exp(100, "_weight_adjust_topo_gbest_100", "gbest", "nonlinear inertia", 10, 1000)
print(i)
i += 1
run_exp(100, "_topo_gbest_100", "gbest", "constriction", 10, 1000)
print(i)
i += 1

run_exp(100, "_weight_adjust_topo_star_100", "star", "nonlinear inertia", 10, 1000)
print(i)
i += 1
run_exp(100, "_topo_star_100", "star", "constriction", 10, 1000)
print(i)
i += 1

run_exp(100, "_weight_adjust_topo_random_100", "rand", "nonlinear inertia", 10, 1000)
print(i)
i += 1
run_exp(100, "_topo_random_100", "rand", "constriction", 10, 1000)
print(i)
i += 1


run_exp(200, "_weight_adjust_200", "lbest", "nonlinear inertia", 10, 1000)
print(i)
i += 1

run_exp(200, "_weight_adjust_topo_gbest_200", "gbest", "nonlinear inertia", 10, 1000)
print(i)
i += 1
run_exp(200, "_topo_gbest_200", "gbest", "constriction", 10, 1000)
print(i)
i += 1

run_exp(200, "_weight_adjust_topo_star_200", "star", "nonlinear inertia", 10, 1000)
print(i)
i += 1
run_exp(200, "_topo_star_200", "star", "constriction", 10, 1000)
print(i)
i += 1

run_exp(200, "_weight_adjust_topo_random_200", "rand", "nonlinear inertia", 10, 1000)
print(i)
i += 1
run_exp(200, "_topo_random_200", "rand", "constriction", 10, 1000)
print(i)
i += 1


run_exp(50, "_weight_adjust_topo_gbest", "gbest", "nonlinear inertia", 10, 1000)
print(i)
i += 1
run_exp(50, "_weight_adjust_topo_star", "star", "nonlinear inertia", 10, 1000)
print(i)
i += 1
run_exp(50, "_weight_adjust_topo_random", "rand", "nonlinear inertia", 10, 1000)
print(i)



