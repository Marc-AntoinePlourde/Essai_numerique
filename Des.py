import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch

n = 100

h_1 = np.linspace(-0.5, 0.5, n)
h_2 = np.linspace(-0.5, -0.49, n)
r, d = 2, 0.3
#Composante contour de gauche
theta = np.linspace(0, np.pi, n)
theta_1, h_1 = np.meshgrid(theta, h_1)
x_1 = (r*np.cos(theta_1))
y_1 = (d+r*np.sin(theta_1))
#Composante contour de gauche
theta_2 = np.linspace(np.pi, 2*np.pi, n)
x_2 = (r*np.cos(theta_2))
y_2 = (-d+r*np.sin(theta_2))
z_1 = h_1
#Composante dessus dessous de gauche
r_1 = np.linspace(0, 2, n)
b_1 = np.outer(r_1, np.cos(theta_2))
c_1 = np.outer(r_1, np.sin(theta_2))
z_2= np.zeros((n, n))

#Composante dessus dessous de droite
b_2 = np.outer(r_1, np.cos(theta))
c_2 = np.outer(r_1, np.sin(theta))

fig = plt.figure()
#Fixer la limites des axes
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-3,3)
ax.set_ylim(-3,3)
ax.set_zlim(-3,3)
#Contour plot
ax.plot_wireframe(x_1, y_1, z_1, rstride = 5, cstride = 5, color = 'k', edgecolors = 'k', alpha = 0.25)
ax.plot_wireframe(x_2, y_2, z_1, rstride = 5, cstride = 5, color = 'k', edgecolors = 'k', alpha = 0.25)
#top et bot gauche
ax.plot_wireframe(b_1, c_1-d, z_2-0.5, rstride = 5, cstride = 5, color = 'k', edgecolors = 'k', alpha = 0.25)
ax.plot_wireframe(b_1, c_1-d, z_2+0.5, rstride = 5, cstride = 5, color = 'k', edgecolors = 'k', alpha = 0.25)
#top et bot droite
ax.plot_wireframe(b_2, c_2+d, z_2-0.5, rstride = 5, cstride = 5, color = 'k', edgecolors = 'k', alpha = 0.25)
ax.plot_wireframe(b_2, c_2+d, z_2+0.5, rstride = 5, cstride = 5, color = 'k', edgecolors = 'k', alpha = 0.25)
# faire des fléche de champ
# Make the grid
x_fleche, y_fleche, z_fleche = np.meshgrid(np.arange(-2.5, 2.5, 0.8), np.arange(-2.5, 2.5, 0.8), np.arange(-0.5, 0))

# Make the direction data for the arrows
u = 0
v = 0
w = 6

ax.quiver(x_fleche, y_fleche, z_fleche, u, v, w, length=1, normalize=True, color='r')
#Nommer les axes
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')


plt.show()