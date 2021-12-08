"""
Assignment 4 Question 1
Completed by: Anmoldeep Singh 180030002
Dated: 21-05-2021
"""

# Importing relevant libraries and modules
import numpy as np
import matplotlib.pyplot as plt
import math as m


# defining the TDMA Algorithm
def tdma(num_terms, low, up, dia, func, sol):
    d1 = np.zeros(num_terms, float)
    rhs1 = np.zeros(num_terms, float)
    d1[0] = dia[0]
    rhs1[0] = func[0]
    for j in range(1, len(dia)):
        d1[j] = dia[j] - (low[j] * up[j - 1] / d1[j - 1])
        rhs1[j] = func[j] - (low[j] * rhs1[j - 1] / d1[j - 1])
    sol[num_terms - 1] = rhs1[num_terms - 1] / d1[num_terms - 1]
    for p in range(len(dia) - 2, -1, -1):
        sol[p] = (rhs1[p] - up[p] * sol[p + 1]) / d1[p]


# Defining required quantities
num = 25
L = 1.0
rho = 1.0
gamma = 0.01
u = 3.0
del_x = L / num
del_t = 0.1
phi_init = 0.5
phi_A = 0.0
phi_B = 1.0
F = rho * u
D = gamma / del_x
P = F / D
mu = (rho * del_x) / del_t  # ap_node
tol = 1.0e-3  # Tolerance

# Defining required arrays
num_sol_old = np.zeros(num)
main_dia = np.zeros(num)
lower_dia = np.zeros(num)
upper_dia = np.zeros(num)
rhs = np.zeros(num)
x_mesh = np.zeros(num)

# Initialising the property phi and generating grid points
for i in range(0, num):
    num_sol_old[i] = phi_init
for i in range(0, num):
    x_mesh[i] = (i + 0.5) * del_x


def CDS():
    # Defining required quantities
    alpha = D - (F / 2)  # ae
    beta = D + (F / 2)  # aw
    delta = alpha + beta + mu  # ap
    converged = False
    num_sol_cds = np.zeros(num)

    # Calculating the numerical solution
    while not converged:
        # Inserting values in diagonals
        for j in range(1, num - 1):
            lower_dia[j] = -beta
            main_dia[j] = delta
            upper_dia[j] = -alpha
            rhs[j] = mu * num_sol_old[j]

        # Defining the boundary conditions
        def apply_bc_dia():
            lower_dia[0] = 0.0
            upper_dia[0] = -alpha
            main_dia[0] = delta + D
            lower_dia[num - 1] = -beta
            upper_dia[num - 1] = 0.0
            main_dia[num - 1] = delta + D
            rhs[0] = mu * num_sol_old[0] + (beta + D) * phi_A
            rhs[num - 1] = mu * num_sol_old[num - 1] + (alpha + D) * phi_B

        # Applying the boundary conditions
        apply_bc_dia()

        # Solving the system
        tdma(num, lower_dia, upper_dia, main_dia, rhs, num_sol_cds)

        # Checking for convergence
        for k in range(0, num):
            t = abs(num_sol_old[k] - num_sol_cds[k])
            if t < tol:
                converged = True

        # Updating the solution
        for j in range(0, num):
            num_sol_old[j] = num_sol_cds[j]

    return num_sol_cds


def Exac():
    # Defining required quantities
    D = gamma / L
    P = F / D
    c = (phi_B - phi_A) / (m.exp(P) - 1)

    # Defining required arrays
    num_sol_exac = np.zeros(num)

    # Calculating the exact solution
    for l in range(0, num):
        num_sol_exac[l] = phi_A + c * (m.exp((P * x_mesh[l]) / L) - 1)

    return num_sol_exac


def Upwind():
    # Defining required quantities
    alpha = D  # ae
    beta = D + F  # aw
    delta = alpha + beta + mu  # ap
    converged = False

    # Defining required arrays
    num_sol_upwind = np.zeros(num)

    # Calculating the numerical solution
    while not converged:
        # Inserting values in diagonals
        for k in range(1, num - 1):
            lower_dia[k] = -beta
            main_dia[k] = delta
            upper_dia[k] = -alpha
            rhs[k] = mu * num_sol_old[k]

        # Defining the boundary conditions
        def apply_bc_dia():
            lower_dia[0] = 0.0
            upper_dia[0] = -alpha
            main_dia[0] = delta + D
            lower_dia[num - 1] = -beta
            upper_dia[num - 1] = 0.0
            main_dia[num - 1] = delta + D
            rhs[0] = mu * num_sol_old[0] + (beta + D) * phi_A
            rhs[num - 1] = mu * num_sol_old[num - 1] + (alpha + D) * phi_B

        # Applying the boundary conditions
        apply_bc_dia()

        # Solving the system
        tdma(num, lower_dia, upper_dia, main_dia, rhs, num_sol_upwind)

        # Checking for convergence
        for k in range(0, num):
            t = abs(num_sol_old[k] - num_sol_upwind[k])
            if t < tol:
                converged = True

        # Updating the solution
        for k in range(0, num):
            num_sol_old[k] = num_sol_upwind[k]

    return num_sol_upwind


def Hybrid():
    # Defining required quantities
    converged = False

    # Defining required arrays
    num_sol_hybrid = np.zeros(num)

    if P < -2:  # Corresponds to Upwind scheme
        alpha = -F  # ae
        beta = 0.0  # aw
    elif -2 <= P <= 2:  # Corresponds to CDS scheme
        alpha = D - F / 2  # ae
        beta = D + F / 2  # aw
    else:   # Corresponds to Upwind scheme
        alpha = 0.0  # ae
        beta = F  # aw

    delta = alpha + beta + mu  # ap

    # Calculating the numerical solution
    while not converged:
        # Inserting values in diagonals
        for p in range(1, num - 1):
            lower_dia[p] = -beta
            main_dia[p] = delta
            upper_dia[p] = -alpha
            rhs[p] = mu * num_sol_old[p]

        # Defining the boundary conditions
        def apply_bc_dia():
            lower_dia[0] = 0.0
            upper_dia[0] = -alpha
            main_dia[0] = delta + D
            lower_dia[num - 1] = -beta
            upper_dia[num - 1] = 0.0
            main_dia[num - 1] = delta + D
            rhs[0] = mu * num_sol_old[0] + (beta + D) * phi_A
            rhs[num - 1] = mu * num_sol_old[num - 1] + (alpha + D) * phi_B

        # Applying the boundary conditions
        apply_bc_dia()

        # Solving the system
        tdma(num, lower_dia, upper_dia, main_dia, rhs, num_sol_hybrid)

        # Checking for convergence
        for k in range(0, num):
            t = abs(num_sol_old[k] - num_sol_hybrid[k])
            if t < tol:
                converged = True

        # Updating the solution
        for p in range(0, num):
            num_sol_old[p] = num_sol_hybrid[p]

    return num_sol_hybrid


# Executing all the schemes
num_sol_Exac = Exac()
num_sol_CDS = CDS()
num_sol_Upwind = Upwind()
num_sol_Hybrid = Hybrid()

# Extending arrays to the boundaries for plotting
num_sol_CDS = np.insert(num_sol_CDS, 0, phi_A)
num_sol_CDS = np.insert(num_sol_CDS, len(num_sol_CDS), phi_B)
num_sol_Upwind = np.insert(num_sol_Upwind, 0, phi_A)
num_sol_Upwind = np.insert(num_sol_Upwind, len(num_sol_Upwind), phi_B)
num_sol_Hybrid = np.insert(num_sol_Hybrid, 0, phi_A)
num_sol_Hybrid = np.insert(num_sol_Hybrid, len(num_sol_Hybrid), phi_B)
num_sol_Exac = np.insert(num_sol_Exac, 0, phi_A)
num_sol_Exac = np.insert(num_sol_Exac, len(num_sol_Exac), phi_B)
x_mesh = np.insert(x_mesh, 0, 0.0)
x_mesh = np.insert(x_mesh, len(x_mesh), L)

# Plotting the final numerically calculated solution
plt.plot(x_mesh, num_sol_CDS, '-bo', label='CDS Solution')
plt.plot(x_mesh, num_sol_Upwind, '-gd', label='Upwind Solution')
plt.plot(x_mesh, num_sol_Hybrid, '-r|', label='Hybrid Solution')
plt.plot(x_mesh, num_sol_Exac, '-c*', label='Exact Solution')
plt.title(chr(961)+' = {0}, '.format(rho)+chr(964)+' = {0}, u = {1}, '.format(gamma, u)+chr(948)+'x = {0}'.format(del_x), fontsize=20)
plt.xlabel('$x$', fontsize=25)
plt.ylabel(chr(966), fontsize=25)
plt.legend(fontsize='large', shadow=True, loc='upper left')
plt.grid()
plt.show()
