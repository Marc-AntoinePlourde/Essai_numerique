import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3


# ATTENTION:
# Ce programme a été très peu testé et est susceptible de comporter de nombreuses erreurs. Utilisez-le à vos risques et périls.

pi = 3.1415926539793238462
t = 0
B_0 = 0.0006
m_0 = 9.10938356*10**(-31)
c = 299792458
q = -1.60217662*10**(-19)
r = 0.01
v_init = np.array([0, 100000, 0])
v = v_init
V = np.linalg.norm(v)
theta = 0
z = 0
f = 0
phi = 0.1
position_de = 0.005
posinit = np.array([0, -0.025, 0])
pos = posinit
iterations = 100000
cadrage = 0.3
r_init = -m_0 * np.linalg.norm(v_init) / (q * B_0)
cadrage_centre = 0 # 0 ou 1
liste = []
delta_t1 = 0.0000006
delta_t2 = delta_t1 / 1000
E_des = np.array([0, 0, 0])
E_entre = np.array([15010, 0, 0])
sauceur_de_premiere = 0
nom_de_fichier = f"dt_{delta_t2}_it_{iterations}"
nom_de_fichier = "fuck_you"
delta = delta_t1
compteur_de_tours = 0
liste_périodes = []
t_1 = 0
def gamma(v):
    """
    Sert à calculer le facteur gamma pour une vitesse donnée.

    param v: vecteur vitesse

    returns:
    facteur gamma pour cette vitesse
    """
    V = np.linalg.norm(v)
    k = 1 / np.sqrt(1 - V ** 2 / c ** 2)
    return k


def champ_electrique(E):
    global t
    #signe du champ
    return np.sign(pos[1]) * E
    # return E * np.sign(np.sin(np.pi * t / 2.4869879291185388e-08))


def champ_magnetique():
    global pos
    # distance à l'orgine
    r = np.linalg.norm([abs(pos[0]) - position_de, pos[1]], pos[2])
    # champ qui dépend du facteur gamma
    if pos[0] >= 2:
        return [0, 0, 0]
    return np.array([0, 0, B_0 / np.sqrt(1 - (r * q * B_0 / (m_0 * c))**2)])


def position():
    global compteur_de_tours
    global t_1
    global E
    global pos
    global t
    global v
    global V
    # masse relativiste
    mgam = gamma(v) * m_0
    if abs(pos[0]) <= position_de:
        # quand la particule se trouve entre les dés
        delta = delta_t1 / np.sqrt(V)
        t += delta
        E = champ_electrique(E_entre)
        a_E =  q * E / (mgam)
        a_B = E_des
    elif (pos + v * delta_t2)[0] * np.sign(pos[0]) < position_de:
        # intermédiaire entre les dés et l'entre-dés
        compteur_de_tours += 0.5
        liste_périodes.append(t - t_1)
        delta = abs((abs(pos[0]) - abs(position_de)) / v[0])
        t_1 = t
        t += delta
        B = champ_magnetique()
        a_B = q * (np.cross(v, B)) / (mgam)
        a_E = E_des
    else:
        # quand la particule se trouve dans les dés
        delta = delta_t2
        t += delta
        # calcul champ électrique
        E = np.array([0,0,0])
        # calcul champ magnétique
        B = champ_magnetique()
        a_B = q * (np.cross(v, B)) / (mgam)
        a_E = E_des
    # calcul module de vitesse
    V = np.linalg.norm(v)
    # calcul vitesse due au champ magnétique
    v_B = v + delta * a_B
    # calcul nouvelle vitesse
    v = (v_B * V / np.linalg.norm(v_B)) + a_E * delta
    # calcul position
    pos = pos + v * delta
    return pos


fig = plt.figure()
# ax = p3.Axes3D(fig)
ax = plt.axes(projection='3d')

nb = 1000




def update(num, data, line):
    if num < 1000:
        n = 0
    else:
        n = num - 1000
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])



# on calcule toutes les positions
for i in range(iterations):
    liste.append(list(position()))



data = np.array(liste).T
line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])


ax.set_xlim3d([posinit[0] * cadrage_centre - r_init - cadrage, posinit[0] * cadrage_centre + r_init + cadrage])
ax.set_xlabel('X')

ax.set_ylim3d([posinit[1] * cadrage_centre - r_init - cadrage, posinit[1] * cadrage_centre + r_init + cadrage])
ax.set_ylabel('Y')

ax.set_zlim3d([posinit[2] * cadrage_centre - r_init - cadrage, posinit[2] * cadrage_centre + r_init + cadrage])
ax.set_zlabel('Z')

ani = animation.FuncAnimation(fig, update, iterations, fargs=(data, line), interval=1, blit=False)
#ani.save(f'{nom_de_fichier}.gif', writer='imagemagick')
#ani.save(f'{nom_de_fichier}.mp4', writer='imagemagick')
# print(liste)
plt.show()

# Données pertinentes à la simulation
print(f"len(liste) = {len(liste)}")
print(f"v = {v}")
V = np.linalg.norm(v)
print(f"demi-période moyenne = {sum(liste_périodes) / len(liste_périodes)}")
print(f"demi-période théorique = {m_0 * np.pi / (q * B_0)}")
print(f"V = {V}")
Beta = V / c
print(f"Beta = {Beta}")
print(f"gamma = {1/np.sqrt(1-Beta**2)}")
print(f"pos = {pos}")
print(f"R final = {np.linalg.norm(pos)}")
print(f"t = {t}")