import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3


pi = 3.1415926539793238462 # pi
c = 299792458 #Vitesse de la lumière en m/s
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


# temps initial
t = 0
# champ magnétique (en Tesla)
m_0 = 1.6726219*10**(-27) #masse au repos utilisée pour le proton
# m_0 = 9.1093837015 * 10**(-31)
q = 1.60217662 * 10**(-19) #charge
v_desiree= 0.70 * c #Vitesse désirée à la fin de l'accélération
r = 2 # rayon des dés
B_0 = (- m_0 * v_desiree) / (q * r) #Champ initial
E_des = np.array([0, 0, 0]) # champ électrique dans les dés (nul)
E_entre = np.array([400000, 0, 0]) # champ entre les dés
v_init = np.array([500000, 0, 0]) # vecteur de vitesse initiale
v = v_init # vitesse
V = np.linalg.norm(v) # grandeur de la vitesse
position_de = 0.05 # Position des dés par rapport à l'axes des x
# rayon de Larmor initial pour calculer la position initiale
r_init = m_0 * np.sqrt(V**2 + 2 * abs(q * np.linalg.norm(E_entre) * position_de / m_0)) / (abs(q) * abs(B_0))
posinit = np.array([0.0000001, - r_init, 0]) # position initiale
pos = posinit # position dans le cyclotron
iterations = 4000000 # nombre d'itérations
liste = [] # liste dans laquelle seront placées toutes les positions
delta_t = 0.00000000006 # pas de temps entre et dans les dés en secondes
delta = delta_t # pas de temps
compteur_de_tours = 0 # un simple nombre qui compte le nombre de tours
liste_periodes = [] # liste dans laquelle seront placées toutes les périodes
t_1 = 0 # temps écoulé depuis que la particule est entrée dans ou sortie des dés
theta = 0


def champ_electrique():
    global t
    global pos
    global v
    # champ qui dépend de la position
    if abs(pos[0]) <= position_de:
        E = E_entre
    else:
        E = E_des
    #signe du champ
    return np.sign(v[0]) * E * np.sign(q)


def champ_magnetique():
    global pos
    global v
    # distance à l'orgine
    r = np.linalg.norm([abs(pos[0]) - position_de, pos[1]])
    # pour déterminer lorsque la particule sort du cyclotron
    if abs(np.linalg.norm(pos)) >= 2 or abs(pos[0]) <= position_de:
        return [0, 0, 0]
    # champ qui dépend du facteur gamma calculé selon la position dans le champ magnétique
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
    if np.sign(abs(pos[0]) - position_de - 0.0001) != np.sign(abs((pos + v * delta_t)[0]) - position_de - 0.0001) and abs(pos[1]) < r:
        # intermédiaire entre les dés et l'entre-dés
        delta = abs((abs(pos[0]) - abs(position_de - 0.0001)) / v[0])

    else:
        # quand la particule se trouve entre ou dans les dés
        delta = delta_t
    t += delta
    # calcul champ électrique et magnétique
    E = champ_electrique()
    B = champ_magnetique()
    # calcul accélération causée par chaque type d'accélération
    a_E = q * E / (mgam)
    a_B = q * np.cross(v, B) / mgam
    # calcul module de vitesse
    V = np.linalg.norm(v)
    # vitesse causée par le champ magnétique
    v_B = v + delta * a_B
    # calcul vitesse finale
    v = (v_B * V / np.linalg.norm(v_B)) + a_E * delta
    #v = v + a_B * delta + a_E * delta
    # calcul position finale
    anc_pos = pos
    pos = pos + v * delta
    if np.sign(pos[0]) != np.sign(anc_pos[0]):
        compteur_de_tours += 0.5
        if pos[0] > 0:
            liste_periodes.append(t - t_1)
            t_1 = t
    return pos


fig = plt.figure()
ax = plt.axes(projection='3d')

nb = 1000




# on calcule toutes les positions
for i in range(iterations):
    liste.append(list(position()))



data = np.array(liste).T
line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])

line.set_data(data[:2, :])
line.set_3d_properties(data[2, :])



ax.set_xlim3d([0.1-r, 0.1+r])
ax.set_xlabel('X')

ax.set_ylim3d([0.1-r, 0.1+r])
ax.set_ylabel('Y')

ax.set_zlim3d([0.1-r,0.1+r])
ax.set_zlabel('Z')

# Données pertinentes à la simulation
print(liste_periodes)
if m_0 == 9.10938356*10**(-31):
    print("particule: électron\n")
elif m_0 == 1.6726219*10**(-27):
    print("particule: proton\n")
else:
    print("particule: inconnu\n")
print(f"iterations = {iterations}")
print(f"delta_t = {delta_t}")
print(f"E = {E_entre}")
print(f"position_de = {position_de}")
print(f"v_désirée = {v_desiree}")
print(f"v_init = {v_init}")
print(f"len(liste) = {len(liste)}")
print(f"v = {v}")
print(f"")
V = np.linalg.norm(v)
print(f"période moyenne = {4 * sum(liste_periodes) / len(liste_periodes)}")
print(f"période théorique = {4 * np.pi * m_0 / (q * B_0)}")
print(f"V = {V}")
Beta = V / c
print(f"Beta = {Beta}")
gamma = 1 / np.sqrt(1-Beta**2)
print(f"gamma = {gamma}")

print(f"pos = {pos}")
print(f"R final = {np.linalg.norm(pos)}")
print(f"t = {t} seconde")
print(f"nbr de tour: {compteur_de_tours} tours")
print(f"Énergie cinétique = {((gamma - 1) * m_0 * c**2) * 6.242 * 10**18} eV")
print(f"Énergie totale = {(gamma * m_0 * c**2) * 6.242 * 10**18 / 1000000} MeV")
print(f"énergie de masse = {m_0 * c**2 * 6.242 * 10**18 / 1000000}")
print(f"Champ magnétique = {B_0} Tesla")


# ani = animation.FuncAnimation(fig, update, iterations, fargs=(data, line), interval=1, blit=False)
#ani.save(f'{nom_de_fichier}.gif', writer='imagemagick')
#ani.save(f'{nom_de_fichier}.mp4', writer='imagemagick')
# print(liste)
plt.show()

