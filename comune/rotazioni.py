import numpy as np

# Strumenti sulle rotazioni comuni ai tre controllori. Convenzioni del repo: z verso l'alto, angoli ZYX
# eta = (phi, theta, psi) con R = Rz(psi) Ry(theta) Rx(phi), omega in terna body.


def angolo_in_pi(a):
    # Riporta un angolo in [-pi, pi]
    return (a + np.pi) % (2*np.pi) - np.pi


def hat(w):
    # Matrice antisimmetrica: hat(w) @ x = w x x
    return np.array([[0, -w[2], w[1]],
                     [w[2], 0, -w[0]],
                     [-w[1], w[0], 0]])


def vee(S):
    # Inversa di hat: dalla matrice antisimmetrica al vettore
    return np.array([S[2, 1], S[0, 2], S[1, 0]])


def matrice_R(eta):
    # Rotazione da terna body a mondo con gli angoli ZYX
    phi, theta, psi = eta
    sf, cf = np.sin(phi), np.cos(phi)
    st, ct = np.sin(theta), np.cos(theta)
    sp, cp = np.sin(psi), np.cos(psi)
    return np.array([[cp*ct, cp*st*sf - sp*cf, cp*st*cf + sp*sf],
                     [sp*ct, sp*st*sf + cp*cf, sp*st*cf - cp*sf],
                     [-st,   ct*sf,            ct*cf           ]])


def angoli_da_R(R):
    # Angoli ZYX dalla matrice, come in Drone.aggiorna_stato
    theta = -np.arcsin(np.clip(R[2, 0], -1, 1))
    phi = np.arctan2(R[2, 1], R[2, 2])
    psi = np.arctan2(R[1, 0], R[0, 0])
    return np.array([phi, theta, psi])


def matrice_W(eta):
    # Cinematica degli angoli ZYX: eta_dot = W(eta) @ omega, con omega in terna body (singolare a theta = +-90°)
    phi, theta, _ = eta
    sf, cf = np.sin(phi), np.cos(phi)
    tt, ct = np.tan(theta), np.cos(theta)
    return np.array([[1, sf*tt, cf*tt],
                     [0, cf,    -sf  ],
                     [0, sf/ct, cf/ct]])
