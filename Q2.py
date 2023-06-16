import numpy as np
from numpy import pi, log, sqrt, array
import matplotlib.pyplot as plt

# Question 2

# Constants
e0 = 8.854187 * (10 ** (-12))
V_up = 1 / 2
V_down = -1 / 2
R = 1
d = 0.025
V1 = 1
D_list = [1 / 2, 1 / 5, 1 / 10, 1 / 50, 1 / 100]
Q_analytical = 8 * e0 * V1 * R


def distance_2D(x_1, y_1, x_2, y_2):
    return ((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2) ** 0.5


def distance_3D(x_1, y_1, z_1, x_2, y_2, z_2):
    return ((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2) ** 0.5


def l_mn_2D(d, x_n, y_n, x_m, y_m):
    return (1 / (4 * pi * e0)) * (d ** 2 / distance_2D(x_n, y_n, x_m, y_m))


def l_mn_3D(d, x_n, y_n, z_n, x_m, y_m, z_m):
    return (1 / (4 * pi * e0)) * (d ** 2 / distance_3D(x_n, y_n, z_n, x_m, y_m, z_m))


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


def l_matrix_AB(R, d, D):
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
                            l_m += [
                                l_mn_3D(d, x_mat[i_n, j_n], y_mat[i_n, j_n], 0, x_mat[i_m, j_m], y_mat[i_m, j_m], D)]
                l_mat += [l_m]
    return array(l_mat)


def capacity(R, d, D, V_up, V_down):
    l_AB = l_matrix_AB(R, d, D)
    l_BA = l_AB.T
    l_AA = l_matrix(R, d)
    l_BB = l_AA
    up_block = np.concatenate((l_AA, l_AB))
    down_block = np.concatenate((l_BA, l_BB))
    l = np.concatenate((up_block, down_block), axis=1)
    V_up_vector = np.ones(len(l_AA[0])) * V_up
    V_down_vector = np.ones(len(l_BB[0])) * V_down
    V = np.concatenate((V_up_vector, V_down_vector))
    sigma = np.linalg.solve(l, V)
    Q_1 = sum(sigma[:len(l_AA)]) * d ** 2
    Capacity = Q_1 / (V_up - V_down)
    return Capacity


CAP = []
for D in D_list:
    CAP.append(capacity(R, d, D, V_up, V_down))

# Calculate theoretical values
C_theo = []
for D in D_list:
    C_theo.append(e0 * pi / D)

# In order not to run the program every time to check the results, we wrote them here:
Answers = [9.450696681991661e-11, 1.831604572670503e-10, 3.289957270906356e-10, 1.5648845868736763e-09, 4.1931163711068e-09,
5.5632497665420493e-11, 1.3908124416355123e-10, 2.7816248832710246e-10, 1.3908124416355123e-09, 2.7816248832710245e-09]

plt.plot(D_list, CAP)
plt.plot(D_list, C_theo)
plt.title("C(D)")
plt.xlabel("D")
plt.ylabel("C [F]")
plt.legend(["C Numerical", "C Theoretical"])
plt.show()

print(Answers)

# D section , Calculate total charge on plates

def total_charge(R, d, D, V_up, V_down):
    l_AB = l_matrix_AB(R, d, D)
    l_BA = l_AB.T
    l_AA = l_matrix(R, d)
    l_BB = l_AA
    up_block = np.concatenate((l_AA, l_AB))
    down_block = np.concatenate((l_BA, l_BB))
    l = np.concatenate((up_block, down_block), axis=1)
    V_up_vector = np.ones(len(l_AA[0])) * V_up
    V_down_vector = np.ones(len(l_BB[0])) * V_down
    V = np.concatenate((V_up_vector, V_down_vector))
    sigma = np.linalg.solve(l, V)
    Q_total = sum(sigma) * d ** 2
    return Q_total

#V up = 1 , V down = 0
print(total_charge(R, d, R / 2, 1, 0))
