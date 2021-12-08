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
num = 200
L = 1.0
rho = 1.0
gamma = 0.01
u = 3.0
del_x = L/num
phi_A = 0.0
phi_B = 1.0
F = rho*u
D = gamma/L
P = F/D
tol = 1.0e-3
# c = (phi_B - phi_A)/(m.exp(P) - 1)
# Defining required arrays
x_mesh = np.zeros(num)
num_sol_new = np.zeros(num)

# Calculating grid points and setting initial values
for i in range(0, num):
    x_mesh[i] = (i + 0.5) * del_x

for l in range(0, num):
    num_sol_new[l] = phi_A + ((phi_B - phi_A)*(m.exp((P*x_mesh[l])/L) - 1))/(m.exp(P) - 1)

# Extending arrays to boundaries for plotting
for i in range(0, 1):
    num_sol_new = np.insert(num_sol_new, 0, phi_A)
    num_sol_new = np.insert(num_sol_new, len(num_sol_new), phi_B)
    x_mesh = np.insert(x_mesh, 0, 0.0)
    x_mesh = np.insert(x_mesh, len(x_mesh), L)
    # x_mesh1 = np.insert(x_mesh1, 0, 0.0)
    # x_mesh1 = np.insert(x_mesh1, len(x_mesh1), L)

# Plotting the final numerically calculated solution
plt.plot(x_mesh, num_sol_new, 'r', label='Numerical Solution @ t=0.1hr')
#
plt.title('Computational solution and Exact solution for %d grid points using TDMA and Fully Implicit Method' % num)
plt.xlabel('$x$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
# plt.legend(fontsize='small', shadow=True, loc='lower right')
plt.grid()
plt.show()
