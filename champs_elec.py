import numpy as np


def champ_electrique(E, temps):
    #signe du champ
    return E*np.sign(np.sin((2*np.pi*m_0)*temps/(q*B_0)))