"""
A run file that sets up Marta's continuous model simulation from the specified config and executes the simulation.
"""

# TODO: Add graph functionality


from sim.continuous_model import matrix_exponential_dist, interaction_matrix_model
import numpy as np
from sim.inputs import get_sim_configuration, get_matrix_function_from_config

cf = get_sim_configuration("continuous")

if cf.subtype == "lattice":
    mat_size = int(cf.no_possible_phenotypes + 1)
    const_affinity_matrix = cf.binding_affinity * np.ones((mat_size, mat_size))
    affinity_matrix = const_affinity_matrix
else:
    get_affinity_matrix = get_matrix_function_from_config(cf.affinity_matrix_config)
    affinity_matrix = 1 * get_affinity_matrix(None)
    print(affinity_matrix)
    # This is deleting rows. Could we do something different?
    real_aff_matrix = np.zeros((21, 21))
    real_aff_matrix[: affinity_matrix.shape[0], : affinity_matrix.shape[1]] = (
        affinity_matrix
    )
    affinity_matrix = real_aff_matrix
    # for i in range(9):
    #    affinity_matrix= np.delete(affinity_matrix,-1, axis=0)


nC, nT, u_vector, time_vector = interaction_matrix_model(
    L=1,
    num_subintervals_lattice=int(cf.no_possible_phenotypes),
    time_step=cf.time_step,
    final_time=cf.final_time,
    A=cf.A,
    a=cf.a,
    theta_c=cf.tumour_selectivity,
    theta_t=cf.CTL_selectivity,
    eta=cf.affinity_range,
    alpha_c=cf.tumour_natural_prolif_rate,
    mu_c=cf.tumour_natural_death_rate,
    zeta_c=cf.tumour_interaction_induced_rate,
    alpha_t=cf.CTL_natural_prolif_rate,
    mu_t=cf.CTL_natural_death_rate,
    zeta_t=cf.CTL_interaction_induced_rate,
    gamma_matrix=affinity_matrix,
    lambda_c=cf.tumour_phenotypic_variation_probability,
    print_results=False,
)

"""
nC, nT, u_vector, time_vector = interaction_matrix_model(
    L=1,
    num_subintervals_lattice=100,
    time_step=0.05,
    final_time=100,
    a=1,
    A=5,
    theta_c=1.8,
    theta_t=1.8,
    eta=2,
    alpha_c=1.5,
    mu_c=5 * 10 ** (-6),
    zeta_c=5 * 10 ** (-6),
    alpha_t=0.05,
    mu_t=5 * 10 ** (-6),
    zeta_t=3 * 10 ** (-5),
    gamma_matrix=const_affinity_matrix,
    lambda_c=0.01,
    print_results=False,
)
"""
