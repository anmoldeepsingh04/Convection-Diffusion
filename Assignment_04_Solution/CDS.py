"""
Question2.1 of Assignment3
Completed by: Anmoldeep Singh 180030002
Dated: 30-04-2021
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
num = 100
L = 1.0
rho = 1.0
gamma = 1.0
u = 1.0
del_x = L/num
del_t = 0.1
phi_init = 0.5
phi_A = 0.0
phi_B = 1.0
F = rho*u
D = gamma/del_x
alpha = D - (F/2)  # ae
beta = D + (F/2)  # aw
mu = (rho*del_x)/del_t  # ap_node
delta = alpha + beta + mu  # ap
converged = False
tol = 1.0e-3
# Defining required arrays
x_mesh = np.zeros(num)
num_sol_old = np.zeros(num)
num_sol_new = np.zeros(num)
main_dia = np.zeros(num)
lower_dia = np.zeros(num)
upper_dia = np.zeros(num)
rhs = np.zeros(num)

# Calculating grid points and setting initial values
for i in range(0, num):
    x_mesh[i] = (i + 0.5) * del_x
    num_sol_old[i] = phi_init

# Calculating the numerical solution
while not converged:
    # Inserting values in diagonals
    for i in range(1, num - 1):
        lower_dia[i] = -beta
        main_dia[i] = delta
        upper_dia[i] = -alpha
        rhs[i] = mu*num_sol_old[i]

    # Defining the boundary conditions
    def apply_bc_dia():
        lower_dia[0] = 0.0
        upper_dia[0] = -alpha
        main_dia[0] = delta + D
        lower_dia[num - 1] = -beta
        upper_dia[num - 1] = 0.0
        main_dia[num - 1] = delta + D
        rhs[0] = mu*num_sol_old[0] + (beta+D)*phi_A
        rhs[num - 1] = mu*num_sol_old[num - 1] + (alpha+D)*phi_B


    # Applying the boundary conditions
    apply_bc_dia()

    # Solving the system
    tdma(num, lower_dia, upper_dia, main_dia, rhs, num_sol_new)

    for k in range(0, num):
        t = abs(num_sol_old[i]-num_sol_new[i])
        if t < tol:
            converged = True

    for i in range(0, num):
        num_sol_old[i] = num_sol_new[i]


# # Calculating the exact solution
# for a in range(0, 1):
#     num1 = 1000
#     del_x1 = 1.0 / num1
#     x_mesh1 = np.zeros(num1)
#     exact_sol_at_1 = np.zeros(num1)
#     exact_sol_at_2 = np.zeros(num1)
#     exact_sol_at_3 = np.zeros(num1)
#     exact_sol_at_4 = np.zeros(num1)
#     exact_sol_at_5 = np.zeros(num1)
#     for i in range(0, num1):
#         x_mesh1[i] = (i + 0.5) * del_x1
#
#     for i in range(0, num1):
#         sum_val = 0.0
#         t = T * (1 / 5)
#         for j in range(1, 10):
#             sum_val += m.exp(-1 * ((j * m.pi) ** 2) * t) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
#         exact_sol_at_1[i] = T_Boundary + 2 * (T_Boundary_0 - T_Boundary) * sum_val
#     for i in range(0, num1):
#         sum_val = 0.0
#         t = T * (2 / 5)
#         for j in range(1, 10):
#             sum_val += m.exp(-1 * ((j * m.pi) ** 2) * t) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
#         exact_sol_at_2[i] = T_Boundary + 2 * (T_Boundary_0 - T_Boundary) * sum_val
#     for i in range(0, num1):
#         sum_val = 0.0
#         t = T * (3 / 5)
#         for j in range(1, 10):
#             sum_val += m.exp(-1 * ((j * m.pi) ** 2) * t) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
#         exact_sol_at_3[i] = T_Boundary + 2 * (T_Boundary_0 - T_Boundary) * sum_val
#     for i in range(0, num1):
#         sum_val = 0.0
#         t = T * (4 / 5)
#         for j in range(1, 10):
#             sum_val += m.exp(-1 * ((j * m.pi) ** 2) * t) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
#         exact_sol_at_4[i] = T_Boundary + 2 * (T_Boundary_0 - T_Boundary) * sum_val
#     for i in range(0, num1):
#         sum_val = 0.0
#         t = T * (5 / 5)
#         for j in range(1, 10):
#             sum_val += m.exp(-1 * ((j * m.pi) ** 2) * t) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
#         exact_sol_at_5[i] = T_Boundary + 2 * (T_Boundary_0 - T_Boundary) * sum_val

# Extending arrays to boundaries for plotting
for i in range(0, 1):
    num_sol_new = np.insert(num_sol_new, 0, phi_A)
    num_sol_new = np.insert(num_sol_new, len(num_sol_new), phi_B)
    x_mesh = np.insert(x_mesh, 0, 0.0)
    x_mesh = np.insert(x_mesh, len(x_mesh), L)
    # x_mesh1 = np.insert(x_mesh1, 0, 0.0)
    # x_mesh1 = np.insert(x_mesh1, len(x_mesh1), L)

# Plotting the final numerically calculated solution
plt.plot(x_mesh, num_sol_new, 'b', label='Numerical Solution @ t=0.1hr')
#
plt.title('Computational solution and Exact solution for %d grid points using TDMA and Fully Implicit Method' % num)
plt.xlabel('$x$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
# plt.legend(fontsize='small', shadow=True, loc='lower right')
plt.grid()
plt.show()

print(x_mesh)