"""
A modification of Marta's simulation, to work seamlessly with my own simulation framework.
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os


def create_u_vector(L, num_subintervals_lattice):
    """
    Creates lattice of points for simulation
    Args
    L: float, (positive) size of interval determined by [-L,L]
    num_points_lattice: int, number of subintervals to divide [-L,L] into
    Returns
    u_vector: array of points in lattice given the number of subintervals interested in
    """
    # Calculate setpsize given interval [-L,L] and number of subintervals we want
    spatial_step = 2 * L / num_subintervals_lattice

    # Create array with the number of points we will have in the lattice
    u_vector = np.zeros(
        num_subintervals_lattice + 1
    )  # Add 1 for right endpoint of lattice

    # Each element in u should be separated by the spatial stepsize
    for i in range(len(u_vector)):
        u_vector[i] = i * spatial_step

    # Subtract L so that the first value of u is -L and the final value is L
    u_vector = u_vector - L

    return u_vector


def initial_conditions(a, A, lattice_points):
    """
    Creates initial densities of cancer cells and T-cells at each of the points in the lattice
    Args:
    a, A: float, parameter values for initial conditions
    lattice_points: array of points in lattice
    returns:
    nC_0: array, initial condition of cancer cell densities
    nT_0: array, initial condition of T-cell densities
    """
    # Initial cancer cells
    nC_0 = (10**4) * (1 + a * np.cos(A * lattice_points))

    # Initial T-cells
    nT_0 = (10**4) * (2 + a * np.cos(A * lattice_points))

    return nC_0, nT_0


def g(x_i, y, chi):
    """
    Defines function g used in intra-population competition and binding terms of model
    Args
    x_i: float, point in lattice
    y: array of all points in lattice
    chi: float, distance determining neighbourhood of x_i that will affect evolution of densities at x_i
    Returns
    g_vals: array of values g takes at each point in the lattice
    """
    # Find difference between point x_i in lattice and each element y_i in y of lattice
    d = abs(x_i - y)

    # Consider the points y_i withing chi of x_i
    close_y_index = np.where(d <= chi)[0]

    # Find size of interval of points within chi of x_i in [-L,L]
    l_max = min(x_i + chi, y[np.max(close_y_index)])
    l_min = max(x_i - chi, y[np.min(close_y_index)])
    l = l_max - l_min

    # Create vector of values of g at each y_i
    g_vals = np.zeros(len(y))

    # If l>0 define the value of g for points withing chi of x_i as 1/l
    if l != 0:
        for index in close_y_index:
            g_vals[index] = 1 / l

    return g_vals


def calculate_competition_vector(lattice_points, selectivity_range, n_vec):
    """
    Calculates using the function g (and depending on the input data) the vectors K_C^h, K_T^h (intra-population competition)
    Args
    lattice_points: array containing all point in lattice
    selectivity_range: float, distance determining neighbourhood of competition
    n_vec: array of densities at each point in lattice of population of cancer cells/T-cells
    Returns
    competition_vector: array containing for each point in the lattice the sum of g applied to all points, multiplied by the corresponding densities and the spatial step
    """
    # Calculate spatial-step
    spatial_step = lattice_points[1] - lattice_points[0]

    # Create vector with zeros of same size as number of lattice points
    competition_vector = np.zeros(len(lattice_points))

    # Iterate through each point in the lattice
    for i in range(len(lattice_points)):
        # Apply the function g given the selectivity range
        g_vec = g(lattice_points[i], lattice_points, selectivity_range)
        # Multiply each element in g_vec by the corresponding element in n_vec and the spatial step, then sum over all elements
        competition_vector[i] = np.sum(spatial_step * np.multiply(g_vec, n_vec))

    return competition_vector


def second_derivative_discretisation(lattice_points):
    """
    Creates matrix for the second central difference used to discretise the second derivative of the cancer cell density (excluding endpoints)
    Args
    lattice_points: array of points in lattice
    Returns
    A: 2-dimensional array
    """
    # Get size of lattice
    n = len(lattice_points)

    # Create central diagonal array
    diagonal = -2 * np.ones(n)
    # Set first and last element to 0
    diagonal[0] = 0
    diagonal[-1] = 0
    # Create upper diagonal array
    upper_diagonal = 1 * np.ones(n - 1)
    # Set first element to 0
    upper_diagonal[0] = 0
    # Create lower diagonal array
    lower_diagonal = 1 * np.ones(n - 1)
    # Set last element to 0
    lower_diagonal[-1] = 0

    # Create Sparse matrix
    diagonals = [lower_diagonal, diagonal, upper_diagonal]
    offset = [-1, 0, 1]
    A = sp.sparse.diags(diagonals, offset)

    return A


def calculate_I(gamma_matrix, lattice_points, selectivity_range, n_vec, spatial_step):
    """
    Calculates the binding by multiplying the interaction matrix by the densities of the population being considered
    Args
    lattice_points: array containing all point in lattice
    selectivity_range: float, distance determining neighbourhood of competition/binding
    n_vec: array of densities at each point in lattice of population of cancer cells/T-cells
    Returns
    binding_vector: array containing for each point in the lattice the sum of g applied to all points, multiplied by the corresponding densities, the interaction matrix and the spatial step
    """
    # Create vector with zeros of size of lattice
    binding_vector = np.zeros(len(lattice_points))

    # Iterate through each point in lattice
    for i in range(len(lattice_points)):
        # Apply the function g given the selectivity range
        g_vec = g(lattice_points[i], lattice_points, selectivity_range)
        # Multiply each element in g_vec by the corresponding element in n_vec, the ith row of the interaction matrix and the spatial step, then sum over all elements
        binding_vector[i] = np.sum(
            spatial_step * np.multiply.reduce([g_vec, n_vec, gamma_matrix[i, :]])
        )

    return binding_vector


def calculate_Rc_interaction_matrix(alpha, mu, K, zeta, Ic):
    """
    Calculates growth rate of cancer cells
    Args
    alpha: float, division rate
    mu: float, intra-population competition death rate
    K: array containing intra-population competition of each cancer cell
    zeta: float, binding death rate
    Ic: array containing binding of each cancer cell
    Returns
    R: array containing the overall growth rate of each cancer cell
    """
    R = alpha - mu * K - zeta * Ic
    return R


def calculate_Rt_interaction_matrix(alpha, mu, K, zeta, It):
    """
    Calculates growth rate of T-cells
    Args
    alpha: float, division rate
    mu: float, intra-population competition death rate
    K: array containing intra-population competition of each T-cell
    zeta: float, binding birth rate
    It: array containing binding of each T-cell
    Returns
    R: array containing the overall growth rate of each T-cell
    """
    R = alpha - mu * K + zeta * It
    return R


def update_parameters(
    nC_k,
    nT_k,
    time_step,
    spatial_step,
    R_c,
    R_t,
    beta_c,
    central_difference_matrix,
):
    """
    Updates cancer cell and T-cell densities at next time-step
    Args
    nC_k: array containing cancer cell densities at each lattice point
    nT_k: array containing T-cell densities at each lattice point
    time_step: float, distance between two consecutive time points at which the model is evaluated
    spatial_step: float, distance between lattice points
    R_c: array containing growth of the cancer cell density at each lattice-point
    R_t: array containing growth of the cancer cell density at each lattice-point
    beta_c: float, diffusion coefficient
    central_difference_matrix: 2-dimensional array, matrix used to discretise the second derivative of the cancer cell densities
    Returns
    nC_k_plus_one: array containing cancer cell densities at each lattice point at the next time-step
    nT_k_plus_one: array containing T-cell densities at each lattice point at the next time-step
    """
    # Cancer cells
    nC_k_plus_one_half = np.multiply(
        nC_k,
        np.divide(
            1 + time_step * np.maximum(0, R_c),
            1 + time_step * abs(np.minimum(0, R_c)),
        ),
    )
    nC_k_plus_one = np.zeros(len(nC_k))
    nC_k_plus_one[1:-1] = (
        nC_k_plus_one_half[1:-1]
        + beta_c
        * (time_step / spatial_step**2)
        * central_difference_matrix.dot(nC_k_plus_one_half)[1:-1]
    )

    # When beta_c=0 we do not implement homogeneous Neumann Boundary Conditions
    if beta_c > 0:
        # Neumann Boundary Conditions
        nC_k_plus_one[0] = nC_k_plus_one[1]
        nC_k_plus_one[-1] = nC_k_plus_one[-2]

    # T-cells
    nT_k_plus_one = np.multiply(
        nT_k,
        np.divide(
            1 + time_step * np.maximum(0, R_t),
            1 + time_step * abs(np.minimum(0, R_t)),
        ),
    )

    return nC_k_plus_one, nT_k_plus_one


def plot_total_cells_evolution(time_vector, total_cancer_cells, total_t_cells):
    """
    Args
    time_vector: array containing each time-step at which the model is evaluated
    total_cancer_cells: array containing the total cancer cells at which model is evaluated
    total_t_cells: array containing the total T-cells at which model is evaluated
    """
    # Plot results of total cells over time
    plt.rcParams[
        "savefig.directory"
    ] = "./outputs"  # os.chdir(os.path.dirname(dir_name))
    plt.title("Evolution of Cancer cells and T-cells over time", fontsize=20)
    plt.plot(time_vector[:-1], total_cancer_cells, label="Cancer cells", color="r")
    plt.plot(time_vector[:-1], total_t_cells, label="T-cells", color="b")
    plt.xlabel("Time", fontsize=18)
    plt.ylabel("Total number of cells", fontsize=18)
    plt.legend(fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.show()


def interaction_matrix_model(
    L,
    num_subintervals_lattice,
    time_step,
    final_time,
    a,
    A,
    theta_c,
    theta_t,
    eta,
    alpha_c,
    mu_c,
    zeta_c,
    alpha_t,
    mu_t,
    zeta_t,
    gamma_matrix,
    lambda_c,
    plot_cell_evolution=True,
    print_results=True,
):
    """
    Simulates evolution of cancer cells and T-cells at each point in [-L,L] being considered over time
    Args
    L: float, (positive) size of interval determined by [-L,L]
    num_subintervals_lattice: int, number of subintervals to divide [-L,L] into
    time_step: float, distance between two consecutive time points at which the model is evaluated
    final_time: float, final time at which model is simulated
    a, A: float, parameter values for initial conditions
    theta_c: float, intra-population competition range of cancer cells
    theta_t: float, intra-population competition range of T-cells
    eta: float, binding distance between cancer cells and T-cells (set to eta=2L unless implementing original version of model)
    alpha_c: float, cancer cell division rate
    mu_c: float, cancer cell intra-population competition death rate
    zeta_c: float, cancer cell binding death rate
    alpha_t: float, T-cell division rate
    mu_t: float, T-cell intra-population competition death rate
    zeta_t: float, T-cell binding birth rate
    gamma_matrix: 2-dimensional array, (i,j)th element contains binding affinity between ith cancer cell and jth T-cell
    lambda_c: float, probability of mutation of cancer cells (used to calculate beta_c)
    plot_cell_evolution: boolean, used to plot total cancer cells and T-cells over time
    print_results: boolean, used to print out total cancer cells and T-cells at each time-step
    Returns
    nC_matrix: list of arrays, each array contains the cancer cell densities across the lattice at a time-step
    nT_matrix: list of arrays, each array contains the T-cell densities across the lattice at a time-step
    u_vector: array containing points in lattice
    tine_vector: array containing all the time points at which the densities are evaluated
    """
    # Calculate spatial step
    spatial_step = 2 * L / num_subintervals_lattice

    # Calculate beta_c
    beta_c = lambda_c * (spatial_step**2) / (2 * time_step)

    # Create empty lists to store cell density at each point over time
    nC_matrix = []
    nT_matrix = []

    # Define lattice points
    u_vector = create_u_vector(L, num_subintervals_lattice)
    v_vector = u_vector.copy()

    # Define initial conditions
    nC_k, nT_k = initial_conditions(a, A, u_vector)
    nC_matrix.append(nC_k)
    nT_matrix.append(nT_k)

    # Define matrix for spatial discretisation
    matrix = second_derivative_discretisation(u_vector)

    # Define total cell lists for each time
    total_cancer_cells = []
    total_t_cells = []

    time = 0
    time_vector = [0]
    # Iterate over each time-step
    while time <= final_time:
        if print_results:
            print("time =", round(time, 2))

        # Count total cells of each type across lattice
        total_cancer_cells.append(np.sum(nC_k) * spatial_step)
        total_t_cells.append(np.sum(nT_k) * spatial_step)

        if time == 0:
            print("Total cancer cells:", total_cancer_cells[-1])
            print("Total T-cells:", total_t_cells[-1])

        # Calculate Kc, Kt, Jc, Jt with the corresponding vectors nC_k and nT_k
        Kc = calculate_competition_vector(u_vector, theta_c, nC_k)
        Kt = calculate_competition_vector(v_vector, theta_t, nT_k)
        gamma_Jc = calculate_I(gamma_matrix, u_vector, eta, nT_k, spatial_step)
        gamma_Jt = calculate_I(
            gamma_matrix.T, v_vector, eta, nC_k, spatial_step
        )  # Input transpose of gamma_matrix so that rows show the binding affinity of each T-cell

        # Use Kc, Kt, Jc, Jt and the matrix gamma to calculate R_c and R_t
        R_c = calculate_Rc_interaction_matrix(alpha_c, mu_c, Kc, zeta_c, gamma_Jc)
        R_t = calculate_Rt_interaction_matrix(alpha_t, mu_t, Kt, zeta_t, gamma_Jt)

        # Update vectors
        nC_k, nT_k = update_parameters(
            nC_k, nT_k, time_step, spatial_step, R_c, R_t, beta_c, matrix
        )
        nC_matrix.append(nC_k)
        nT_matrix.append(nT_k)

        # Update time
        time += time_step
        time_vector.append(time)

        # Print results for total cells at each timestep
        if print_results:
            print("Total cancer cells:", total_cancer_cells[-1])
            print("Total T-cells:", total_t_cells[-1])

    # Plot results for total cells over time
    if plot_cell_evolution:
        plot_total_cells_evolution(time_vector, total_cancer_cells, total_t_cells)

    return nC_matrix, nT_matrix, u_vector, time_vector


def matrix_exponential_dist(num_subintervals_lattice, delta):
    """
    Using an exponential distribution, randomly generates the entries of the interaction matrix
    Args
    num_subintervals_lattice: int, number of subintervals to divide [-L,L] into
    delta: float, parameter of exponential distribution
    """
    matrix = np.random.exponential(
        scale=delta,
        size=(num_subintervals_lattice + 1, num_subintervals_lattice + 1),
    )
    return matrix


run = False

if __name__ == "__main__" and run:
    # Randomly generate an interaction matrix
    gamma_mat = matrix_exponential_dist(num_subintervals_lattice=100, delta=1)

    # Implement simulation
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
        gamma_matrix=gamma_mat,
        lambda_c=0.01,
    )
