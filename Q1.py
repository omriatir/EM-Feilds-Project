import numpy as np
from numpy import pi, log, sqrt, array
import matplotlib.pyplot as plt

# Question 1

# Constants
e0 = 8.854187 * (10 ** (-12))
V = 1
R = 1
d_list = [0.25, 0.15, 0.12, 0.1, 0.075, 0.05, 0.025, 0.02]

Q_analytical = 8 * e0 * V * R

def distance_2D(x_1, y_1, x_2, y_2):
    return ((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2) ** 0.5


def l_mn_2D(d, x_n, y_n, x_m, y_m):
    return (1 / (4 * pi * e0)) * (d ** 2 / distance_2D(x_n, y_n, x_m, y_m))


def l_nn(d):
    return d * log((sqrt(2) + 1) / (2 ** 0.5 - 1)) / (2 * pi * e0)


def x_matrix(R, d):
    a = int(2 * R // d)
    return array([[d * i + d / 2 for i in range(a)] for j in range(a)])


def y_matrix(R, d):
    a = int(2 * R // d)
    return array([[d * j + d / 2 for i in range(a)] for j in range(a)])


def l_matrix(R, d):
    a = int(2 * R // d)
    x_mat = x_matrix(R, d)
    y_mat = y_matrix(R, d)
    l_mat = []
    for i_n in range(a):
        for j_n in range(a):
            if distance_2D(x_mat[i_n, j_n], y_mat[i_n, j_n], R, R) < R:
                l_m = []
                for i_m in range(a):
                    for j_m in range(a):
                        if distance_2D(x_mat[i_m, j_m], y_mat[i_m, j_m], R, R) < R:
                            if (i_n, j_n) == (i_m, j_m):
                                l_m += [l_nn(d)]
                            else:
                                l_m += [l_mn_2D(d, x_mat[i_n, j_n], y_mat[i_n, j_n], x_mat[i_m, j_m], y_mat[i_m, j_m])]
                l_mat += [l_m]

    return array(l_mat)

l_list = []
sigma_list = []
Q_list = []

#Calculate l_matrix
for d in d_list:
    l_list.append(l_matrix(R,d))

#Solve Linear System
for i in range(len(l_list)):
    sigma_list.append(np.linalg.solve(l_list[i],np.ones(len(l_list[i]))))
    Q_list.append(sum(sigma_list[i])*d_list[i]**2)

print(Q_list)

# In order not to run the program every time to check the results, we wrote them here:

ANSWERS = [7.045845327598558e-11, 7.038071292396779e-11, 6.959561412821214e-11, 6.928680008859111e-11,
           7.005294791722533e-11, 7.043496196223577e-11, 7.058598054386094e-11, 7.069292361889084e-11]

plt.plot(d_list,ANSWERS)
plt.plot(d_list,np.ones(len(d_list))*Q_analytical)
plt.title("Q(d)")
plt.xlabel("d")
plt.ylabel("Q [c]")
plt.legend(["Q Numerical","Q Analytical"])
plt.show()



