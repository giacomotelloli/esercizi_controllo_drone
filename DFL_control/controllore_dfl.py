import numpy as np

from comune.rotazioni import angolo_in_pi, matrice_R, matrice_W

E3 = np.array([0, 0, 1.0])

# S @ w = w x e3 = (w_y, -w_x, 0): serve nelle derivate di b3 = R e3
S_E3 = np.array([[0, 1.0, 0],
                 [-1.0, 0, 0],
                 [0, 0, 0]])


def derivata_W3(eta, eta_dot):
    # Derivata nel tempo della terza riga di W(eta), quella che dà psi_dot = W3 omega
    phi, theta, _ = eta
    dphi, dtheta, _ = eta_dot
    sf, cf = np.sin(phi), np.cos(phi)
    st, ct = np.sin(theta), np.cos(theta)
    return np.array([0,
                     (cf*dphi*ct + sf*st*dtheta) / ct**2,
                     (-sf*dphi*ct + cf*st*dtheta) / ct**2])


class ControlloreDFL:
    """
    Linearizzazione in retroazione dinamica (slide 03 del corso, modello completo), con z verso l'alto.

    Uscite (x, y, z, psi). Il compensatore zeta = (T, T') rende T'' un ingresso; derivando quattro volte
    la posizione e due volte lo yaw compare l'ingresso u~ = (T'', omega'):
        [p''''; psi''] = l(xi, zeta) + A(xi, zeta) [T''; omega']        tau = J omega' + omega x J omega
    Si sceglie [T''; omega'] = A^-1 (v - l), quindi p'''' = v_p e psi'' = v_psi: catene di integratori
    chiuse con i coefficienti c3..c0 (posizione) e c1, c0 (yaw).

    Come nel corso, il modello è scritto sugli angoli di Eulero: R si ricostruisce da eta, e A è singolare
    per T = 0 e per theta = +-90° (dove W non esiste).
    """

    def __init__(self, massa:float, J, guadagni:dict, dt:float, g:float = 9.81, limiti:dict = None):

        # Parametri del modello usato dal controllore (come nel BS si usa solo la diagonale di J)
        self.m = massa
        self.g = g
        self.J = np.diag(np.diag(np.array(J, dtype=float)))

        # Coefficienti dei polinomi d'errore: righe x, y, z con colonne (c3, c2, c1, c0); yaw (c1, c0)
        self.C = np.array([guadagni['x'], guadagni['y'], guadagni['z']], dtype=float)
        self.c_psi = np.array(guadagni['psi'], dtype=float)

        # Protezione dalla singolarità a spinta nulla (planning §5.7): il compensatore non scende sotto T_min.
        # Cambia il metodo e va dichiarata; di default lo stesso 10% di mg del limite su F_z del BS
        limiti = limiti or {}
        self.T_min = limiti.get('T_min', 0.1*massa*g)

        self.dt = dt
        self.reset()

        return

    def reset(self, T0:float = None):
        # Compensatore all'hover: zeta(0) = (mg, 0), oppure la spinta iniziale data
        self.zeta = np.array([self.m*self.g if T0 is None else T0, 0.0])

        # Grandezze interne dell'ultimo passo, lette dal main per il log
        self.e_p = np.zeros(3)          # errore di posizione p_d - p [m]
        self.e_psi = 0.0                # errore di yaw [rad]
        self.v_virt = np.zeros(4)       # ingressi virtuali (v_p, v_psi)
        self.det_A = 0.0                # determinante di A: vicino a zero ci si avvicina alla singolarità
        self.eta_d = np.full(3, np.nan) # il DFL non ha un assetto desiderato: solo lo yaw (per i grafici)

        return

    def calcola(self, stato:dict, rif:dict):
        # Stato misurato + riferimento -> spinta totale T [N] e coppie tau [N*m] in terna body
        # stato: p, v, eta, omega (R non serve: il modello è sugli angoli di Eulero)
        # rif: p, v, a, j, s (fino allo snap), psi, dpsi, ddpsi
        p, v, eta, omega = stato['p'], stato['v'], stato['eta'], stato['omega']
        m, g = self.m, self.g
        T, dT = self.zeta

        R = matrice_R(eta)
        b3 = R[:, 2]
        W = matrice_W(eta)
        eta_dot = W @ omega

        # Derivate dell'uscita dal modello: accelerazione e jerk usano la spinta del compensatore, non una misura
        w_e3 = S_E3 @ omega                             # omega x e3
        acc = T/m * b3 - g*E3                           # p''
        jerk = dT/m * b3 + T/m * R @ w_e3               # p'''

        # Ingressi virtuali: errore p_d - p e sue derivate, chiusi con i polinomi scelti
        self.e_p = rif['p'] - p
        errori = [rif['j'] - jerk, rif['a'] - acc, rif['v'] - v, self.e_p]     # e''', e'', e', e
        v_p = rif['s'] + sum(self.C[:, k]*errori[k] for k in range(4))

        self.e_psi = angolo_in_pi(rif['psi'] - eta[2])
        v_psi = rif['ddpsi'] + self.c_psi[0]*(rif['dpsi'] - eta_dot[2]) + self.c_psi[1]*self.e_psi
        self.v_virt = np.concatenate([v_p, [v_psi]])

        # Modello: p'''' = T''/m b3 + 2 T'/m R (omega x e3) + T/m R (omega x (omega x e3)) + T/m R (omega' x e3)
        #          psi'' = W3' omega + W3 omega'
        l = np.zeros(4)
        l[:3] = 2*dT/m * R @ w_e3 + T/m * R @ np.cross(omega, w_e3)
        l[3] = derivata_W3(eta, eta_dot) @ omega

        A = np.zeros((4, 4))
        A[:3, 0] = b3/m
        A[:3, 1:] = T/m * R @ S_E3
        A[3, 1:] = W[2]
        self.det_A = np.linalg.det(A)

        u = np.linalg.solve(A, self.v_virt - l)
        ddT, domega = u[0], u[1:]
        tau = self.J @ domega + np.cross(omega, self.J @ omega)

        self.eta_d = np.array([np.nan, np.nan, rif['psi']])

        # La spinta comandata è lo stato del compensatore; poi lo si integra (Eulero in avanti) con il nuovo T''
        T_cmd = T
        self.zeta = self.zeta + self.dt*np.array([dT, ddT])
        if self.zeta[0] < self.T_min:
            self.zeta[0] = self.T_min
            self.zeta[1] = max(self.zeta[1], 0.0)   # non continua a scendere contro il limite

        return T_cmd, tau
