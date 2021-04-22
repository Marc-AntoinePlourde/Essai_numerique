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
posinit = np.array([0.0009476050108216742, 0, 0])
pos = posinit
print(pos)
iterations = 100000
cadrage = 0.3
r_init = -m_0 * np.linalg.norm(v_init) / (q * B_0)
cadrage_centre = 0 # 0 ou 1
liste = []
delta_t1 = 0.0000006
delta_t2 = delta_t1 / 1000
E_des = np.array([0, 0, 0])
E_entre = np.array([30000, 0, 0])
sauceur_de_premiere = 0
nom_de_fichier = f"dt_{delta_t2}_it_{iterations}"
nom_de_fichier = "fuck_you"
delta = delta_t1
print(m_0 * np.linalg.norm(v) / (q * B_0))
print(f"r = {m_0 * np.linalg.norm(v) / (q * B_0)}")
compteur_de_tours = 0
theta = pi / 2 - 1.08430149339592
B_0 = B_0 / np.sqrt(1 - np.cos(theta))
print(np.tan(theta))
t_1 = 0
def gamma(v):
    # global sauceur_de_premiere
    V = np.linalg.norm(v)
    # print(f"V = {V}")
    #try:
    k = 1 / np.sqrt(1 - V ** 2 / c ** 2)

    #except RuntimeWarning:
    #    sauceur_de_premiere = 1
    return k


def addition(u_prime, v):
    global pos
    u_primerad = u_prime * (u_prime @ v) / (v @ v)
    tan = np.cross(np.array([0,0,1]), v)
    u_primetan = u_prime * (u_prime @ tan) / (tan @ tan)
    u_rad = (u_primerad+v)/(1 + v @ u_primerad/c**2)
    u_tan = u_primetan / ((1+v @ u_primerad /c**2) * gamma(v))
    return u_rad + u_tan


def transfo(pos, v, t):
    """
    Transformation de Lorentz de la position et du temps dans le référentiel du cyclotron vers celui de la particule.
    Pour avoir la transformation inverse, simplement mettre un signe négatif à v.
    :param pos: np.array de la position (m)
    :param v: np.array de la vitesse de la particule par rapport au référentiel de l'observatoire (m/s)
    :param t: temps en secondes
    :return: liste du temps et de la position dans le référentiel primé
    """
    gam = gamma(v)
    V = np.linalg.norm(v)
    n = (1 / V) * v
    t_prime = gam*(t-V*n @ pos/c**2)
    r_prime = pos + (gam - 1)*(pos@n)*n-gam*t*V*n
    return [t_prime, r_prime]


def invtransfo(r_prime, v, t_prime):
    V = np.linalg.norm(v)
    n = (1 / V) * v
    gam = gamma(v)
    t = gam * (t_prime+r_prime @ n * V/c**2)
    r = r_prime + (gam - 1)*(r_prime @ n) * n + gam*t_prime*V*n
    return [t, r]


def transfelec(E, B, v):
    """

    Applique la transformée de Lorentz au champ électrique et au champ magnétique.

    :param E: np.array du champ électrique
    :param B: np.array du champ magnétique
    :param v: np.array de la vitesse entre les deux référentiels
    :return: liste du champ électrique et du champ magnétique dans le référentiel primé.
    """
    gam = gamma(v)
    n = v / np.linalg.norm(v)
    E_prime = gam * (E + np.cross(v, B)) - (gam - 1) * (E @ n) * n
    B_prime = gam * (B - np.cross(v, E) / c**2) - (gam - 1) * (B @ n) * n
    return [E_prime, B_prime]


def transfacc(a, v, u):
    """
    Applique la transformée de Lorentz de l'accélération.
    :param a: np.array de l'accélération (m/s**2)
    :param v: np.array de la vitesse entre les deux référentiels (m/s)
    :param u: np.array de la vitesse de la particule dans le référentiel observé (m/s)
    :return: l'accélération dans le référentiel primé
    """
    gam = gamma(v)
    gammo = gam * (1 - v @ u / c**2)
    V = np.linalg.norm(v)
    return a / gammo**2 - (gam - 1) * (a @ v) * v / (V**2*gammo**3) + (a @ v) * u * gam/(c**2*gammo**3)


def champ_electrique(E):
    global t
    global f
    #signe du champ
    # if t > 1000 * delta_t2:
        #return E_des
    return np.sign(pos[1]) * E
    # return E * np.sign(np.sin((q * B_0) * t / m_0))


def champ_magnetique():
    global pos
    # distance à l'orgine
    r = np.linalg.norm([abs(pos[0]) - position_de, pos[1]])
    # champ qui dépend du facteur gamma
    #if abs(pos[1] / pos[0]) < 0.5288941338363227:
    if pos[0] >= 2:
        return [0, 0, 0]
    return np.array([0, 0, B_0 / np.sqrt(1 - (r * q * B_0 / (m_0 * c))**2)])


liste_de_pisse = []
def position():
    global liste_de_pisse
    global compteur_de_tours
    global E
    global pos
    global t
    global v
    global V
    global f
    global sauceur_de_premiere
    global t_1
    # masse relativiste
    mgam = gamma(v) * m_0
    if abs(pos[0]) <= position_de:
        delta = delta_t1 / np.sqrt(V)
        t += delta
        E = champ_electrique(E_entre)
        a_E =  q * E / (mgam)
        a_B = E_des
    elif (pos + v * delta_t2)[0] * np.sign(pos[0]) < position_de:
        if compteur_de_tours == 0.5:
           print(f"R = {m_0 * np.linalg.norm(v) / (q * B_0)}")
        liste_de_pisse.append(t - t_1)
        #if len(liste_de_pisse) == 10:
            #print(liste_de_pisse)
            #liste_de_pisse = []
        compteur_de_tours += 0.5
        delta = abs((abs(pos[0]) - abs(position_de)) / v[0])
        t_1 = t
        t += delta
        B = champ_magnetique()
        a_B = q * (np.cross(v, B)) / (mgam)
        a_E = E_des
    else:
        delta = delta_t2
        t += delta
        E = np.array([0,0,0])
        B = champ_magnetique()
        a_B = q * (np.cross(v, B)) / (mgam)
        a_E = E_des
    V = np.linalg.norm(v)
    v_B = v + delta * a_B
    v = (v_B * V / np.linalg.norm(v_B)) + a_E * delta
    pos = pos + v * delta
    # v = v_B * V / np.linalg.norm(v_B) + a_E * delta
    return pos


fig = plt.figure()
# ax = p3.Axes3D(fig)
ax = plt.axes(projection='3d')

nb = 1000



# def animate(i):
#     global pos
#     global vit
#     pos = position()[0]
#     vit = position()[1]
#     return(pos, vit)

def update(num, data, line):
    if num < 1000:
        n = 0
    else:
        n = num - 1000
    line.set_data(data[:2, n:num])
    line.set_3d_properties(data[2, n:num])




for i in range(iterations):
    # print(f"i = {i}")
    # print(f"t = {t}")
    # if E[0] != 0 or sauceur_de_premiere == 1:
        # break
    # if i == 481:
        # v = np.array([30000000000000000000, 98, 0])
    liste.append(list(position()))



# for i in range(int(iterations/10)):
#    liste.append(liste[-1])

data = np.array(liste).T
line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])


ax.set_xlim3d([posinit[0] * cadrage_centre - r_init - cadrage, posinit[0] * cadrage_centre + r_init + cadrage])
ax.set_xlabel('X')

ax.set_ylim3d([posinit[1] * cadrage_centre - r_init - cadrage, posinit[1] * cadrage_centre + r_init + cadrage])
ax.set_ylabel('Y')

ax.set_zlim3d([posinit[2] * cadrage_centre - r_init - cadrage, posinit[2] * cadrage_centre + r_init + cadrage])
ax.set_zlabel('Z')

print("Will tu sôces.")
# print(data)
ani = animation.FuncAnimation(fig, update, iterations, fargs=(data, line), interval=1, blit=False)
print("fuck you")
#ani.save(f'{nom_de_fichier}.gif', writer='imagemagick')
#ani.save(f'{nom_de_fichier}.mp4', writer='imagemagick')
# print(liste)
plt.show()

# new = position()

# pos = new[0]
# vit = new[1]

# liste.append(pos)
print(f"len(liste) = {len(liste)}")
print(f"v = {v}")
V = np.linalg.norm(v)
print(f"moyenne = {sum(liste_de_pisse) / len(liste_de_pisse)}")
print(f"V = {V}")
Beta = V / c
print(f"Beta = {Beta}")
print(f"gamma = {1/np.sqrt(1-Beta**2)}")
print(f"pos = {pos}")
print(f"R = {np.linalg.norm(pos)}")
print(m_0 * V / (q * B_0))
print(f"t = {t}")
print(0.1386204919909253 / 0.6901018201639635)
print(f"nombre de tours: {compteur_de_tours}")