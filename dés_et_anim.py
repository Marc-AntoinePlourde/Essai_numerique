import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3


# ATTENTION:
# Ce programme a été très peu testé et est susceptible de comporter de nombreuses erreurs. Utilisez-le à vos risques et périls.

pi = 3.1415926539793238462
# temps initial
t = 0

# champ magnétique (en Tesla)
#m_0 = 1.6726219*10**(-27) #masse au repos utilisée
m_0 = 9 * 10**(-31)
c = 299792458 #Vitesse de la lumière en m/s
q = -1.60217662 * 10**(-19) #charge
v_desire= 0.70 * c #Vitesse désirée à la fin de l'accélération
r = 2 # rayon des dés
B_0 = (- m_0 * v_desire) / (q * r) #Champ initial
E_des = np.array([0, 0, 0])
E_entre = np.array([15000, 0, 0])
v_init = np.array([100000, 0, 0])
v = v_init
V = np.linalg.norm(v)
position_de = 0.1 # Depuis l'axe des x
r_init = m_0 * np.sqrt(V**2 + 2 * abs(q) * np.linalg.norm(E_entre) * position_de / m_0) / (abs(q) * abs(B_0))
posinit = np.array([0, - r_init / 2, 0])
pos = posinit
iterations = 100000 #nombre d'itération fait
cadrage_centre = 0 # 0 ou 1
liste = []
delta_t1 = 0.0000006
delta_t2 = delta_t1 / 1000

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
    return - np.sign(pos[1]) * E *np.sign(q)
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
        V = np.linalg.norm(v)
        v_B = v + delta * a_B
        v = (v_B * V / np.linalg.norm(v_B)) + a_E * delta
        pos = pos + v * delta
    elif (pos + v * delta_t2)[0] * np.sign(pos[0]) < position_de + 0.0001:
        # intermédiaire entre les dés et l'entre-dés
        compteur_de_tours += 0.5
        liste_périodes.append(t - t_1)
        delta = abs((abs(pos[0]) - abs(position_de)) / v[0])
        t_1 = t
        t += delta
        B = champ_magnetique()
        a_B = q * (np.cross(v, B)) / (mgam)
        a_E = E_des
        V = np.linalg.norm(v)
        signe = np.sign(v[0])
        v_B = v + delta * a_B
        v = (v_B * V / np.linalg.norm(v_B)) + a_E * delta
        pos = pos + v * delta
        v = V * signe * np.array([1, 0, 0])
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
        V = np.linalg.norm(v)
        v_B = v + delta * a_B
        v = (v_B * V / np.linalg.norm(v_B)) + a_E * delta
        pos = pos + v * delta
    # calcul module de vitesse

    # calcul vitesse due au champ magnétique

    # calcul nouvelle vitesse

    # calcul position
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
#Pour les dés
n = 100

h_1 = np.linspace(-0.5, 0.5, n)
h_2 = np.linspace(-0.5, -0.49, n)

#Composante contour de gauche
thetaa = np.linspace(np.pi / 2, 3 * np.pi / 2, n)
theta_1, h_1 = np.meshgrid(thetaa, h_1)
x_1 = r * np.cos(theta_1) - position_de
y_1 = r * np.sin(theta_1)
#Composante contour de droite
theta_2 = np.linspace(-np.pi / 2, np.pi / 2, n)
x_2 = position_de + r * np.cos(theta_2)
y_2 = r * np.sin(theta_2)
z_1 = h_1
#Composante dessus dessous de gauche
r_1 = np.linspace(0, 2, n)
b_1 = np.outer(r_1, np.cos(theta_2))
c_1 = np.outer(r_1, np.sin(theta_2))
z_2= np.zeros((n, n))

#Composante dessus dessous de droite
b_2 = np.outer(r_1, np.cos(thetaa))
c_2 = np.outer(r_1, np.sin(thetaa))
#Contour plot
ax.plot_wireframe(x_1, y_1, z_1, rstride = 5, cstride = 5, color = 'k', edgecolors = 'k', alpha = 0.15)
ax.plot_wireframe(x_2, y_2, z_1, rstride = 5, cstride = 5, color = 'k', edgecolors = 'k', alpha = 0.15)
#top et bot gauche
ax.plot_wireframe(b_1+position_de, c_1, z_2-0.5, rstride = 5, cstride = 5, color = 'k', edgecolors = 'k', alpha = 0.15)
ax.plot_wireframe(b_1+position_de, c_1, z_2+0.5, rstride = 5, cstride = 5, color = 'k', edgecolors = 'k', alpha = 0.15)
#top et bot droite
ax.plot_wireframe(b_2-position_de, c_2, z_2-0.5, rstride = 5, cstride = 5, color = 'k', edgecolors = 'k', alpha = 0.15)
ax.plot_wireframe(b_2-position_de, c_2, z_2+0.5, rstride = 5, cstride = 5, color = 'k', edgecolors = 'k', alpha = 0.15)

ax.set_xlim3d([posinit[0] * cadrage_centre - (r + 0.1), posinit[0] * cadrage_centre + (r+0.1)])
ax.set_xlabel('X')

ax.set_ylim3d([posinit[1] * cadrage_centre - (r+0.1), posinit[1] * cadrage_centre + (r+0.1)])
ax.set_ylabel('Y')

ax.set_zlim3d([posinit[2] * cadrage_centre - (r+0.1), posinit[2] * cadrage_centre + (r+0.1)])
ax.set_zlabel('Z')

ani = animation.FuncAnimation(fig, update, iterations, fargs=(data, line), interval=1, blit=False)
#ani.save(f'{nom_de_fichier}.gif', writer='imagemagick')
#ani.save(f'{nom_de_fichier}.mp4', writer='imagemagick')
# print(liste)
plt.show()

# Données pertinentes à la simulation
print(f"len(liste) = {len(liste)}")
print(f"v = {v}")
print(f"")
V = np.linalg.norm(v)
print(f"demi-période moyenne = {sum(liste_périodes) / len(liste_périodes)}")
print(f"demi-période théorique = {m_0 * np.pi / (q * B_0)}")
print(f"V = {V}")
Beta = V / c
print(f"Beta = {Beta}")
gamma = 1 / np.sqrt(1-Beta**2)
print(f"gamma = {gamma}")

print(f"pos = {pos}")
print(f"R final = {np.linalg.norm(pos)}")
print(f"t = {t} seconde")
print(f"nbr de tour: {compteur_de_tours} tours")
print(f"Énergie cinétique = {((gamma - 1) * m_0 * c**2) / abs(q)} eV")
print(f"Énergie totale = {(m_0 * np.sqrt(V**4 * gamma ** 2  + c**4)) / (1000000 * abs(q))} MeV")
print(f"Champ magnétique = {B_0} Tesla")
