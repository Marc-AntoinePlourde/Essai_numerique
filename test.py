import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3


# ATTENTION:
# Ce programme a été très peu testé et est susceptible de comporter de nombreuses erreurs. Utilisez-le à vos risques et périls.


t = 0
B_0 = 0.0002
m_0 = 9.10938356*10**(-31)
c = 299792458
q = -1.60217662*10**(-19)
r = 0.01
v_init = np.array([0, 100000, 0])
v = v_init
theta = 0
phi = 270
z = 0
position_de = 0.002
posinit = np.array([0.0025, 0, 0])
pos = posinit
print(pos)
iterations = 10000
cadrage = 0.02
r_init = -m_0 * np.linalg.norm(v_init) / (q * B_0)
cadrage_centre = 0 # 0 ou 1
liste = []
delta_t = 0.000000002
E = np.array([0, 0, 0])
sauceur_de_premiere = 0
nom_de_fichier = f"dt_{delta_t}_it_{iterations}"

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




def position():
    global E
    global pos
    global t
    global v
    global sauceur_de_premiere
    t += delta_t
    if abs(pos[0]) < position_de:
        E = np.array([20, 0, 0]) * np.sign(pos[1])
    else:
        E = np.array([0,0,0])
    mgam = gamma(v) * m_0
    # print(f"pos = {pos}")
    # if pos[0] == pos[1]:
    #    sauceur_de_premiere = 1
    # print(f"pos[0] = {pos[0]}")
    # print(f"pos[1] = {pos[1]}")
    # r = np.sqrt(pos[0]**2+pos[1]**2)
    # print(r)
    # print(f"sauce = {1 - (q * B_0 * r / (m_0 * c))**2}")
    B = np.array([0, 0, B_0 * gam])
    # * gamma(v)])
    # E_prime, B_prime = transfelec(E, B, v)
    # print(f"E = {E}")
    # print(f"B = {B}")
    a_E = q * E / (mgam)
    a_B = q * (np.cross(v, B)) / (mgam)
    # print(f"a_E = {a_E}")
    # print(f"a_B = {a_B}")
    # cross = np.cross(v, B)
    # print(f"cross = {np.cross(v, B)}")
    # print(f"norme cross = {np.linalg.norm(cross)}")
    # print(f"dot prod = {v @ cross}")
    # print(f"a_prime = {a_prime}")
    # r_prime, t_prime = transfo(pos, v, delta_t)
    # print(f"r_prime = {r_prime}")
    # print(f"t_prime = {t_prime}")
    V = np.linalg.norm(v)
    # print(f"V = {np.linalg.norm(v)}")
    pos = pos + v * delta_t
    v_B = (v + delta_t * a_B)
    v = v_B * V / np.linalg.norm(v_B) + a_E * delta_t
    # print(f"v_f = {v_f}")
    # r_prime2 = r_prime + t_prime * v_f
    # v = addition(v_f, v)
    # print(f"v = {v}")
    # transfo(r_prime2, -u, t_prime)
    return pos


fig = plt.figure()
# ax = p3.Axes3D(fig)
ax = plt.axes(projection='3d')

nb = 50




# def animate(i):
#     global pos
#     global vit
#     pos = position()[0]
#     vit = position()[1]
#     return(pos, vit)

def update(num, data, line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])




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


# print(data)
ani = animation.FuncAnimation(fig, update, nb, fargs=(data, line), interval=1000 / nb, blit=False)
ani.save(f'{nom_de_fichier}.gif', writer='imagemagick')
ani.save(f'{nom_de_fichier}.mp4', writer='imagemagick')
# print(liste)
plt.show()

# new = position()

# pos = new[0]
# vit = new[1]

# liste.append(pos)
print(f"len(liste) = {len(liste)}")
print(f"v = {v}")
V = np.linalg.norm(v)
print(f"V = {V}")
print(f"pos = {pos}")
print(m_0 * V / (q * B_0))