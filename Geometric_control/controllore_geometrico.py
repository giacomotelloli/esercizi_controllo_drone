import numpy as np

from comune.rotazioni import hat, vee, angoli_da_R

E3 = np.array([0, 0, 1.0])


class ControlloreGeometrico:
    """
    Controllo geometrico su SE(3) (Lee, Leok, McClamroch, CDC 2010), riscritto con z verso l'alto:
    m p'' = -m g e3 + f R e3. L'assetto si controlla direttamente su SO(3), senza angoli di Eulero.

    Guadagni diagonali fisici, uno per asse (con guadagni scalari, come nell'articolo, quota e yaw non
    potrebbero avere poli diversi da x, y e da rollio, beccheggio):
        kx, kv [N/m, N*s/m]           -> poli di posizione  s^2 + (kv/m) s + kx/m
        kR, kOmega [N*m/rad, N*m*s]   -> poli d'assetto     s^2 + (kOmega/J) s + kR/J
    Vanno ricalcolati con m e J della scena: se restano quelli nominali, nello scenario B i poli si spostano.

    Omega_d e Omega_d': tau_filtro=None -> nulli; altrimenti derivata filtrata di R_d, come nel BS
    (stessa scelta per i due, STEPS 0.6).

    direzione_yaw: come si costruisce R_d dallo yaw desiderato
      'zyx' (default): b1d nel piano verticale della direzione di yaw, quindi lo yaw ZYX di R_d è esattamente
             psi_d, come per BS, DFL e la traiettoria (costruzione di Mellinger e Kumar)
      'lee': b1d è la proiezione di (cos psi_d, sin psi_d, 0) sul piano normale a b3d, come nelle slide;
             con il drone inclinato lo yaw ZYX di R_d differisce da psi_d
    limiti: {'T_min': ...} spinta minima, di default 0.1 mg come negli altri due. Con f = F . R e3 < 0
            (drone oltre 90° da b3d) il mixer del Drone porterebbe i motori a zero e toglierebbe anche la coppia.
    """

    def __init__(self, massa:float, J, guadagni:dict, dt:float, g:float = 9.81, tau_filtro = None,
                 direzione_yaw:str = 'zyx', limiti:dict = None):

        # Parametri del modello usato dal controllore (come negli altri due si usa solo la diagonale di J)
        self.m = massa
        self.g = g
        self.J = np.diag(np.diag(np.array(J, dtype=float)))

        self.kx = np.diag(guadagni['kx'])
        self.kv = np.diag(guadagni['kv'])
        self.kR = np.diag(guadagni['kR'])
        self.kOmega = np.diag(guadagni['kOmega'])

        if direzione_yaw not in ('zyx', 'lee'):
            raise ValueError("direzione_yaw: 'zyx' oppure 'lee'")
        self.direzione_yaw = direzione_yaw
        limiti = limiti or {}
        self.T_min = limiti.get('T_min', 0.1*massa*g)

        self.dt = dt
        self.tau_filtro = tau_filtro

        self.reset()

        return

    def reset(self):
        # Memoria per le derivate di R_d
        self.primo_passo = True
        self.Rd_prec = np.eye(3)
        self.dRd = np.zeros((3, 3))
        self.Omega_d = np.zeros(3)
        self.dOmega_d = np.zeros(3)
        self.b2d_prec = np.array([0, 1.0, 0])

        # Grandezze interne dell'ultimo passo, lette dal main per il log
        self.e_p = np.zeros(3)      # errore di posizione p_d - p [m] (stesso segno degli altri controllori)
        self.R_d = np.eye(3)        # assetto desiderato
        self.eta_d = np.zeros(3)    # suoi angoli ZYX, solo per i grafici
        self.e_R = np.zeros(3)      # errore d'assetto su SO(3)
        self.e_Omega = np.zeros(3)  # errore di velocità angolare [rad/s]
        self.Psi = 0.0              # funzione d'errore d'assetto 1/2 tr(I - R_d^T R), in [0, 2]

        return

    def calcola(self, stato:dict, rif:dict):
        # Stato misurato + riferimento -> spinta f [N] e momenti M [N*m] in terna body
        # stato: p, v, R, omega; rif: p, v, a, psi (jerk e snap non servono con la derivata filtrata)
        p, v, R, omega = stato['p'], stato['v'], stato['R'], stato['omega']

        # Posizione: forza desiderata e sua proiezione sull'asse attuale del drone
        e_x = p - rif['p']
        e_v = v - rif['v']
        self.e_p = -e_x
        F = -self.kx @ e_x - self.kv @ e_v + self.m*self.g*E3 + self.m*rif['a']
        f = max(F @ R[:, 2], self.T_min)

        # Assetto desiderato: b3d lungo F, b1d dalla direzione di yaw desiderata
        b3d = F / np.linalg.norm(F)
        sp, cp = np.sin(rif['psi']), np.cos(rif['psi'])
        if self.direzione_yaw == 'lee':
            b2d = np.cross(b3d, np.array([cp, sp, 0]))          # b1d = proiezione di b1c sul piano normale a b3d
        else:
            b1d = np.cross(np.array([-sp, cp, 0]), b3d)         # b1d nel piano verticale di (cos psi, sin psi)
            b2d = np.cross(b3d, b1d)
        n = np.linalg.norm(b2d)
        b2d = b2d/n if n > 1e-6 else self.b2d_prec     # b3d parallelo alla direzione usata: si tiene il precedente
        self.b2d_prec = b2d
        b1d = np.cross(b2d, b3d)
        Rd = np.column_stack([b1d, b2d, b3d])
        self.R_d = Rd
        self.eta_d = angoli_da_R(Rd)

        Omega_d, dOmega_d = self._derivate_assetto(Rd)

        # Errori su SO(3) e momento
        self.e_R = 0.5*vee(Rd.T @ R - R.T @ Rd)
        self.e_Omega = omega - R.T @ Rd @ Omega_d
        self.Psi = 0.5*np.trace(np.eye(3) - Rd.T @ R)

        M = (-self.kR @ self.e_R - self.kOmega @ self.e_Omega
             + np.cross(omega, self.J @ omega)
             - self.J @ (hat(omega) @ R.T @ Rd @ Omega_d - R.T @ Rd @ dOmega_d))

        return f, M

    def _derivate_assetto(self, Rd):
        # Omega_d e Omega_d' da R_d, con lo stesso filtro del BS: s/(tau_f s + 1), Eulero all'indietro
        if self.primo_passo:
            # Nessun valore precedente: si parte con derivate nulle, senza picco iniziale
            self.Rd_prec = Rd.copy()
            self.primo_passo = False

        if self.tau_filtro is None:
            self.Omega_d = np.zeros(3)
            self.dOmega_d = np.zeros(3)
        else:
            a = self.tau_filtro / (self.tau_filtro + self.dt)
            b = 1 / (self.tau_filtro + self.dt)
            self.dRd = a*self.dRd + b*(Rd - self.Rd_prec)
            A = Rd.T @ self.dRd
            Omega_d = vee(0.5*(A - A.T))        # parte antisimmetrica: R_d^T R_d' lo è solo a meno del filtro
            self.dOmega_d = a*self.dOmega_d + b*(Omega_d - self.Omega_d)
            self.Omega_d = Omega_d

        self.Rd_prec = Rd.copy()

        return self.Omega_d, self.dOmega_d
