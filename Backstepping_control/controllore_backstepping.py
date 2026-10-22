import numpy as np


def angolo_in_pi(a):
    # Riporta un angolo in [-pi, pi]
    return (a + np.pi) % (2*np.pi) - np.pi


def matrice_W(eta):
    # Cinematica degli angoli ZYX: eta_dot = W(eta) @ omega, con omega in terna body
    phi, theta, _ = eta
    sf, cf = np.sin(phi), np.cos(phi)
    tt, ct = np.tan(theta), np.cos(theta)
    return np.array([[1, sf*tt, cf*tt],
                     [0, cf,    -sf  ],
                     [0, sf/ct, cf/ct]])


class ControlloreBackstepping:

    def __init__(self, massa:float, J, guadagni:dict, limiti:dict, dt:float, g:float = 9.81, tau_filtro = None):

        # Parametri fisici del modello usato dal controllore (per i test di robustezza possono essere diversi da quelli veri)
        self.m = massa
        self.g = g
        J = np.array(J)
        self.Ix, self.Iy, self.Iz = J[0,0], J[1,1], J[2,2]  # il progetto sul modello semplificato usa solo la diagonale

        # Guadagni (k1, k2) di ogni canale: righe x,y,z e phi,theta,psi
        self.k_pos = np.array([guadagni['x'], guadagni['y'], guadagni['z']])
        self.k_att = np.array([guadagni['phi'], guadagni['theta'], guadagni['psi']])

        # Limiti
        self.angolo_max = limiti['angolo_max']  # inclinazione massima della spinta [rad]
        self.T_max = limiti['T_max']            # spinta massima concessa [N]

        # Derivate di phi_d e theta_d: None -> poste a zero, altrimenti derivata filtrata con costante di tempo tau_filtro [s]
        self.dt = dt
        self.tau_filtro = tau_filtro

        self.reset()

        return

    def reset(self):
        # Memoria per le derivate di (phi_d, theta_d)
        self.primo_passo = True
        self.ang_d_prec = np.zeros(2)
        self.dang_d = np.zeros(2)
        self.ddang_d = np.zeros(2)

        # Grandezze interne dell'ultimo passo, lette dal main per il log
        self.mu = np.zeros(3)       # accelerazioni desiderate [m/sec^2]
        self.eta_d = np.zeros(3)    # angoli desiderati [rad]
        self.e_p = np.zeros(3)      # errore di posizione [m]
        self.e_eta = np.zeros(3)    # errore d'assetto [rad]
        self.V_pos = 0.0            # Lyapunov del sottosistema di posizione
        self.V_att = 0.0            # Lyapunov del sottosistema d'assetto

        return

    def calcola(self, p, v, eta, omega, rif:dict):
        # Stato misurato + riferimento -> spinta totale T [N] e coppie tau [N*m] in terna body
        eta_dot = matrice_W(eta) @ omega

        self.mu = self._posizione(p, v, rif)
        T, phi_d, theta_d = self._inversione(self.mu, eta[2])

        eta_d = np.array([phi_d, theta_d, rif['psi']])
        deta_d, ddeta_d = self._derivate_assetto(eta_d, rif)

        tau = self._assetto(eta, eta_dot, eta_d, deta_d, ddeta_d)

        return T, tau

    def _canale(self, e, de, acc_d, k):
        # Ricetta del backstepping per doppi integratori x1_dot = x2, x2_dot = u (lavora su più canali insieme)
        k1, k2 = k[:,0], k[:,1]
        u = acc_d + (1 + k1*k2)*e + (k1 + k2)*de

        e2 = de + k1*e                      # secondo errore: x2d - x2
        V = 0.5*np.sum(e**2 + e2**2)        # Lyapunov dei canali

        return u, V

    def _posizione(self, p, v, rif):
        # Backstepping su x, y, z con le accelerazioni come ingressi virtuali
        self.e_p = rif['p'] - p
        de_p = rif['v'] - v
        mu, self.V_pos = self._canale(self.e_p, de_p, rif['a'], self.k_pos)

        return mu

    def _inversione(self, mu, psi):
        # Forza desiderata in terna mondo: (T/m) R e3 = (mu_x, mu_y, mu_z + g)
        F = self.m * np.array([mu[0], mu[1], mu[2] + self.g])

        # La spinta non può tirare verso il basso
        F[2] = max(F[2], 0.1*self.m*self.g)

        # Limite d'inclinazione: la componente orizzontale non supera F_z*tan(angolo_max)
        F_oriz = np.linalg.norm(F[:2])
        F_oriz_max = F[2]*np.tan(self.angolo_max)
        if F_oriz > F_oriz_max:
            F[:2] = F[:2] * F_oriz_max / F_oriz

        norma = np.linalg.norm(F)
        T = min(norma, self.T_max)

        # Angoli desiderati: la direzione della spinta, ruotata con lo yaw misurato
        sp, cp = np.sin(psi), np.cos(psi)
        phi_d = np.arcsin(np.clip((F[0]*sp - F[1]*cp)/norma, -1, 1))
        theta_d = np.atan2(F[0]*cp + F[1]*sp, F[2])

        return T, phi_d, theta_d

    def _derivate_assetto(self, eta_d, rif):
        # Derivate prima e seconda di phi_d e theta_d; quelle di psi_d arrivano dalla traiettoria
        ang_d = eta_d[:2]

        if self.primo_passo:
            # Nessun valore precedente: si parte con derivate nulle, senza picco iniziale
            self.ang_d_prec = ang_d.copy()
            self.primo_passo = False

        if self.tau_filtro is None:
            self.dang_d = np.zeros(2)
            self.ddang_d = np.zeros(2)
        else:
            # Derivata "sporca" s/(tau_f s + 1), discretizzata con Eulero all'indietro
            a = self.tau_filtro / (self.tau_filtro + self.dt)
            b = 1 / (self.tau_filtro + self.dt)
            dang_d = a*self.dang_d + b*(ang_d - self.ang_d_prec)
            self.ddang_d = a*self.ddang_d + b*(dang_d - self.dang_d)
            self.dang_d = dang_d

        self.ang_d_prec = ang_d.copy()

        deta_d = np.array([self.dang_d[0], self.dang_d[1], rif['dpsi']])
        ddeta_d = np.array([self.ddang_d[0], self.ddang_d[1], rif['ddpsi']])

        return deta_d, ddeta_d

    def _assetto(self, eta, eta_dot, eta_d, deta_d, ddeta_d):
        # Backstepping su phi, theta, psi sul modello semplificato
        self.eta_d = eta_d
        self.e_eta = eta_d - eta
        self.e_eta[2] = angolo_in_pi(self.e_eta[2])     # niente giri completi sullo yaw
        de_eta = deta_d - eta_dot

        nu, self.V_att = self._canale(self.e_eta, de_eta, ddeta_d, self.k_att)  # accelerazioni angolari desiderate

        # Compensazione dei termini giroscopici del corpo rigido
        dphi, dtheta, dpsi = eta_dot
        tau_phi = self.Ix*nu[0] - (self.Iy - self.Iz)*dtheta*dpsi
        tau_theta = self.Iy*nu[1] - (self.Iz - self.Ix)*dphi*dpsi
        tau_psi = self.Iz*nu[2] - (self.Ix - self.Iy)*dphi*dtheta

        return np.array([tau_phi, tau_theta, tau_psi])
