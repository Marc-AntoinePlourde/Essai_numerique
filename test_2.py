import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.animation import FuncAnimation



# ATTENTION:
# Ce programme a été très peu testé et est susceptible de comporter de nombreuses erreurs. Utilisez-le à vos risques et périls.


t = 0
B_0 = 0.0001
m_0 = 9.10938356*10**(-31)
c = 299792458
q = -1.60217662*10**(-19)
r = 0.1
v = np.array([0, 10000, 0])
theta = 0
phi = 90
z = 0
pos = np.array([r*np.cos(theta), r*np.sin(theta), z])
print(pos)
liste = [pos]
delta_t = 0.000000000002
E = np.array([0, 0, 0])
sauceur_de_premiere = 0


def gamma(v):
    global sauceur_de_premiere
    k = 1
    V = np.linalg.norm(v)
    print(f"V = {V}")
    try:
        k = 1 / np.sqrt(1 - V ** 2 / c ** 2)

    except RuntimeWarning:
        sauceur_de_premiere = 1
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
    if abs(pos[0]) < 0.002:
        E = np.array([-2, 0, 0]) * np.sign(pos[1])
    else:
        E = np.array([0,0,0])
    # print(f"pos = {pos}")
    # if pos[0] == pos[1]:
    #    sauceur_de_premiere = 1
    # print(f"pos[0] = {pos[0]}")
    # print(f"pos[1] = {pos[1]}")
    # r = np.sqrt(pos[0]**2+pos[1]**2)
    # print(r)
    # print(f"sauce = {1 - (q * B_0 * r / (m_0 * c))**2}")
    B = np.array([0, 0, B_0])
    # * gamma(v)])
    # E_prime, B_prime = transfelec(E, B, v)
    # print(f"E = {E}")
    # print(f"B = {B}")
    a = q * (E + np.cross(v, B)) / (m_0 * gamma(v))
    # print(f"a = {a}")
    # print(f"cross = {np.cross(-v, B_prime)}")
    # print(f"a_prime = {a_prime}")
    # r_prime, t_prime = transfo(pos, v, delta_t)
    # print(f"r_prime = {r_prime}")
    # print(f"t_prime = {t_prime}")
    # print(f"V = {np.linalg.norm(v)}")
    pos = pos + v * delta_t
    v = v + delta_t * a
    # print(f"v_f = {v_f}")
    # r_prime2 = r_prime + t_prime * v_f
    # v = addition(v_f, v)
    # print(f"v = {v}")
    # transfo(r_prime2, -u, t_prime)
    return pos

fig = plt.figure()
# ax = p3.Axes3D(fig)
ax = plt.axes(projection='3d')
ln = ax.plot(pos[0], pos[1], pos[2])
nb = 50


def init():
    ax.set_zlim(-6, 6)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim(-a, a)
    ax.set_ylim(-b, b)
    return ln,

def update(frame):
    global ln
    global pos
    global v
    global E
    global t
    #ln.remove()
    t = frame*0.01
    print(f" pos[0] = {pos[0]}")
    position()

    ln = ax.scatter(pos[0], pos[1], pos[2], c='r')
    return ln



ani = FuncAnimation(fig, update, frames=np.arange(frame_depart, frame_fin), interval=20,
                    init_func=init, blit=True)

ani.save(f'graphique_de_pisse.gif', writer='imagemagick')
#plt.show()
