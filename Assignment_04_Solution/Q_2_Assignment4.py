"""
Question2.2 of Assignment4
Completed by: Anmoldeep Singh 180030002
Dated: 23-05-2021
"""

# Importing relevant libraries and modules
import timeit  # To calculate the execution time
import numpy as np
import scipy.integrate as ig
import matplotlib.pyplot as plt

start = timeit.default_timer()  # Starting the timer as code starts from here


# defining the Hybrid Scheme and passing relevant arguments
def Hybrid(peclet, a_w, a_e, a_n, a_s, ap, f_we, f_ns, d_we, d_ns):
    # Setting the coefficients according to the Peclet number
    for m in range(0, len(peclet)):
        if peclet[m] < -2:  # Corresponds to Upwind scheme
            a_e[m] = -f_we[m]  # ae
            a_w[m] = 0.0  # aw
            a_n[m] = -f_ns[m]  # an
            a_s[m] = 0.0  # as
        elif -2 <= peclet[m] <= 2:  # Corresponds to CDS scheme
            a_e[m] = d_we[m] - f_we[m] / 2  # ae
            a_w[m] = d_we[m] + f_we[m] / 2  # aw
            a_n[m] = d_ns[m] - f_ns[m] / 2  # an
            a_s[m] = d_ns[m] + f_ns[m] / 2  # as
        elif peclet[m] > 2:  # Corresponds to Upwind scheme
            a_e[m] = 0.0  # ae
            a_w[m] = f_we[m]  # aw
            a_n[m] = 0.0  # an
            a_s[m] = f_ns[m]  # as

        ap[m] = a_e[m] + a_w[m] + a_n[m] + a_s[m]  # ap


# Defining required quantities
H = 1.0
L = 20.0
del_x = 0.01
del_y = 0.01
num_x = 100
num_y = 100
k = 1.0
c = 100.0
rho_o = 1.0
rho = rho_o * c
T_wall = 100.0
T_entrance = 50.0
converged = False
error = 1.0e-3
v = 0.0
umean = 1.0

# Defining required arrays
x_mesh = np.zeros(num_x)
y_mesh = np.zeros(num_y)
y_mesh_center = np.linspace(-0.5, 0.5, num_y)
num_sol = np.zeros((num_y, num_x))
num_sol_old = np.zeros((num_y, num_x))
F_we = np.zeros(num_y)
F_ns = np.zeros(num_y)
D_we = np.zeros(num_y)
D_ns = np.zeros(num_y)
u = np.zeros(num_y)
Pe = np.zeros(num_y)
aE = np.zeros(num_y)
aW = np.zeros(num_y)
aN = np.zeros(num_y)
aS = np.zeros(num_y)
aP = np.zeros(num_y)
Prod_vel_temp = np.zeros(num_x)
T_Bulk = np.zeros(num_x)
h = np.zeros(num_x)
Nu = np.zeros(num_x)

# Calculating grid points
for i in range(0, num_x):
    x_mesh[i] = i * del_x
for i in range(0, num_y):
    y_mesh[i] = i * del_y

# Calculating the Peclet number and the velocity at a cross section
for i in range(0, num_y):
    u[i] = 1.5 * (1 - 4 * (y_mesh_center[i] ** 2))
    Pe[i] = (rho * u[i] * del_x) / k

# Calculating the Diffusive and Convective coefficients
for i in range(0, num_y):
    F_we[i] = rho * u[i] * del_y
    F_ns[i] = rho * v * del_x
    D_we[i] = (k * del_y) / del_x
    D_ns[i] = (k * del_x) / del_y

# Calling the Hybrid scheme
Hybrid(Pe, aW, aE, aN, aS, aP, F_we, F_ns, D_we, D_ns)

# Defining a counter to count the number of iterations
iter_count = 0

"""Calculating the numerical solution"""

# Defining and applying the boundary conditions
for i in range(0, num_x):
    num_sol[0][i] = 100.0
    num_sol[num_y - 1][i] = 100.0
for k in range(0, num_y):
    num_sol[k][0] = 50.0

while not converged:
    iter_count = iter_count + 1
    print(iter_count)
    # Employing Point by Point Gauss-Seidel method
    for i in range(1, num_x - 1):
        for j in range(1, num_y - 1):
            num_sol[j][i] = (aW[j] * num_sol[j][i - 1] + aN[j] * num_sol[j + 1][i] + aS[j] * num_sol[j - 1][i] + aE[j] *
                             num_sol[j][i + 1]) / aP[j]

    # Checking for convergence
    t = 0.0  # Stores the maximum error in each iteration
    for i in range(1, num_y - 2):
        for j in range(1, num_x - 2):
            temp = abs(num_sol_old[i][j] - num_sol[i][j])
            if temp > t:
                t = temp
    if t < error:
        converged = True

    # Updating the solution
    for i in range(0, num_y - 1):
        for j in range(0, num_x - 1):
            num_sol_old[i][j] = num_sol[i][j]

# Employing the boundary condition on the outflow boundary
for i in range(0, num_y - 1):
    num_sol[i][num_x - 1] = num_sol[i][num_x - 2]

for i in range(0, num_x):
    for j in range(0, num_y):
        Prod_vel_temp[i] = (num_sol[j][i] * u[i])/(H*umean)

# Calculating the bulk mean temperature
for i in range(0, num_x):
    T_Bulk[i] = ig.trapz(y_mesh, Prod_vel_temp)

# Calculating the heat transfer coefficient
for i in range(0, num_x):
    h[i] = (k * (T_wall - num_sol[num_y - 1][i])) / ((T_wall - T_Bulk[i]) * del_y)

# Calculating the Nusselt number
for i in range(0, num_x):
    Nu[i] = (h[i] * 2 * H) / k

# Program execution ends here. Displaying the time of execution and number of iterations
stop = timeit.default_timer()
execution_time = stop - start

print("Program Executed in " + str(execution_time), "seconds")
print("Program Executed took ", iter_count, "iterations")


# Plotting the final exact solution

# Plots the 2D graph of temperature values at all the nodes
c = plt.imshow(num_sol, cmap='hot', vmin=num_sol.min(), vmax=num_sol.max(), extent=[0, 20, 0, 1],
               interpolation='bicubic', origin='lower', aspect='auto')

# # Plots the contour of temperature values in the domain
# c = plt.contourf(x_mesh, y_mesh, num_sol, 25)

plt.colorbar()
plt.title('Numerical Solution after {0} iterations'.format(iter_count), fontweight="bold")
plt.xlabel('$Spatial Dimension$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
plt.get_current_fig_manager().window.state('zoomed')
plt.show()

# Plotting the values of h, Nu
plt.plot(x_mesh, h, label="h")
plt.plot(x_mesh, Nu, label="Nu")
plt.xlabel('$x$', fontsize=15)
plt.ylabel('h and Nu', fontsize=15)
plt.legend(fontsize='large', shadow=True, loc='upper left')
plt.title('Variation of h and Nu')
plt.legend()
plt.show()


