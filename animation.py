import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3

# -----Définition des variables initiales---------



# ---------Modélisation du Réacteur nucléaire-----------
# th = np.linspace(0,2. * np.pi,, stop)
# x = (R_int + R_bob * np.cos())

t = 0


pos = [np.cos(np.pi*t), np.sin(np.pi * t), t]

liste = [pos]
def position():
    global pos
    global t
    t += 0.1
    delta_theta = 0.1
    pos = [np.cos(np.pi*t), np.sin(np.pi*t), t]
    return pos


fig = plt.figure()
ax = p3.Axes3D(fig)

nb = 2000




# def animate(i):
#     global pos
#     global vit
#     pos = position()[0]
#     vit = position()[1]
#     return(pos, vit)

def update(num, data, line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])


for i in range(50000):
    liste.append(position())

data = np.array(liste).T
line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])

ax.set_xlim3d([-15, 15])
ax.set_xlabel('X')

ax.set_ylim3d([-15, 15])
ax.set_ylabel('Y')

ax.set_zlim3d([-15, 15])
ax.set_zlabel('Z')

ani = animation.FuncAnimation(fig, update, nb, fargs=(data, line), interval=1000 / nb, blit=False)
# ani.save('matplot003.gif', writer='imagemagick')
plt.show()

# new = position()

# pos = new[0]
# vit = new[1]

# liste.append(pos)
