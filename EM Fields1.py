import numpy as np
import matplotlib.pyplot as plt

a=3
R=2
L=3
y = np.arange(-3,3.1,0.1)
theta = np.arctan(y/a)
r = a/np.cos(theta)

r_new = R**2/r
x_image = r_new*np.cos(theta)
y_image = r_new*np.sin(theta)
theta1 = np.linspace(0, 2 * np.pi, 100)

x = R * np.cos(theta1)
y = R * np.sin(theta1)

plt.plot(x,y,x_image,y_image)
plt.show()