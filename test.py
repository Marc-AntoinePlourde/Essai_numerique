import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3


# ATTENTION:
# Ce programme a été très peu testé et est susceptible de comporter de nombreuses erreurs. Utilisez-le à vos risques et périls.


t = 0
B_0 = 1
m_0 = 9.10938356*10**(-31)
c = 299792458
q = 1.60217662*10**(-19)
r = 0.1
v = np.array([0, 100, 0])
theta = 0
phi = 90
z = 0
pos = np.array([r*np.cos(theta), r*np.sin(theta), z])
print(pos)
liste = [pos]
delta_t = 0.01
