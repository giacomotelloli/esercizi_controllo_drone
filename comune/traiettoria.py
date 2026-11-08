import numpy as np
from math import factorial


# Metodo di interpolazione -> ordine k della derivata minimizzata: min ∫ ||q^(k)(t)||^2 dt
# La soluzione è un polinomio di grado 2k-1 per tratto, con derivate continue fino alla 2k-2 nei punti interni
METODI = {
    'lineare': 1,   # grado 1: velocità a gradini, accelerazione infinita nei punti (solo per confronto)
    'cubica': 2,    # grado 3: accelerazione continua, jerk a gradini
    'min_jerk': 3,  # grado 5: continua fino allo snap
    'min_snap': 4,  # grado 7: continua fino alla derivata sesta
}


def base(n_coef, r, tau):
    # Derivata r-esima della base 1, tau, tau^2, ... valutata in tau
    riga = np.zeros(n_coef)
    for i in range(r, n_coef):
        riga[i] = factorial(i) / factorial(i-r) * tau**(i-r)
    return riga


class TraiettoriaInterpolata:

    def __init__(self, punti, metodo:str = 'min_snap', tempi = None, v_media:float = None, t_inizio:float = 0.0):
        # punti: N righe (x, y, z) oppure (x, y, z, psi) [m, rad]; senza psi lo yaw desiderato resta 0
        # tempi: istante di passaggio di ogni punto [s], oppure v_media [m/s] per ricavarli dalle distanze

        if metodo not in METODI:
            raise ValueError("Metodo sconosciuto: "+metodo+", scegli tra "+str(list(METODI)))
        if (tempi is None) == (v_media is None):
            raise ValueError("Passa i tempi oppure v_media (uno solo dei due)")

        punti = np.array(punti, dtype=float)
        if len(punti) < 2:
            raise ValueError("Servono almeno due punti (per un hover ripeti lo stesso punto)")
        if punti.shape[1] == 3:
            punti = np.column_stack([punti, np.zeros(len(punti))])
        punti[:,3] = np.unwrap(punti[:,3])  # lo yaw ruota sempre dal verso più corto
        self.punti = punti

        if tempi is not None:
            self.tempi = np.array(tempi, dtype=float)
        else:
            distanze = np.linalg.norm(np.diff(punti[:,:3], axis=0), axis=1)
            self.tempi = t_inizio + np.concatenate([[0], np.cumsum(distanze / v_media)])

        self.durate = np.diff(self.tempi)   # durata di ogni tratto [s]
        if len(self.tempi) != len(punti) or np.any(self.durate <= 0):
            raise ValueError("Serve un tempo per punto, strettamente crescente (con v_media: due punti consecutivi coincidono)")

        self.metodo = metodo
        self.k = METODI[metodo]
        self.t_inizio = self.tempi[0]
        self.t_fine = self.tempi[-1]

        # Coefficienti (tratti, 2k, 4) nella variabile normalizzata tau = (t - t_j)/T_j in [0, 1]
        self.coef = self.calcola_coefficienti()

        return

    def calcola_coefficienti(self):
        # Un unico sistema lineare A c = b per tutti i tratti; le quattro colonne di b sono x, y, z, psi
        k = self.k
        n_coef = 2*k
        n_tratti = len(self.durate)
        T = self.durate

        A = np.zeros((n_tratti*n_coef, n_tratti*n_coef))
        b = np.zeros((n_tratti*n_coef, 4))
        riga = 0

        # 1. Ogni tratto passa per i suoi due punti
        for j in range(n_tratti):
            A[riga, j*n_coef:(j+1)*n_coef] = base(n_coef, 0, 0)
            b[riga] = self.punti[j]
            A[riga+1, j*n_coef:(j+1)*n_coef] = base(n_coef, 0, 1)
            b[riga+1] = self.punti[j+1]
            riga += 2

        # 2. Nei punti interni le derivate da 1 a 2k-2 sono continue
        # d^r/dt^r = (1/T_j^r) d^r/dtau^r: l'uguaglianza è moltiplicata per T_j^r per tenere il sistema ben condizionato
        for j in range(n_tratti-1):
            for r in range(1, 2*k-1):
                A[riga, j*n_coef:(j+1)*n_coef] = base(n_coef, r, 1)
                A[riga, (j+1)*n_coef:(j+2)*n_coef] = -(T[j]/T[j+1])**r * base(n_coef, r, 0)
                riga += 1

        # 3. Partenza e arrivo da fermo: derivate da 1 a k-1 nulle
        for r in range(1, k):
            A[riga, :n_coef] = base(n_coef, r, 0)
            A[riga+1, -n_coef:] = base(n_coef, r, 1)
            riga += 2

        coef = np.linalg.solve(A, b)

        return coef.reshape(n_tratti, n_coef, 4)

    def __call__(self, t:float):
        # Riferimento all'istante t, nello stesso formato del dict rif dei controllori
        if t <= self.t_inizio or t >= self.t_fine:
            # Prima della partenza e dopo l'arrivo: fermo nel primo o nell'ultimo punto
            q = self.punti[0].copy() if t <= self.t_inizio else self.punti[-1].copy()
            d = [q] + [np.zeros(4) for _ in range(4)]
        else:
            j = np.searchsorted(self.tempi, t) - 1     # tratto che contiene t
            tau = (t - self.tempi[j]) / self.durate[j]
            d = [base(2*self.k, r, tau) @ self.coef[j] / self.durate[j]**r for r in range(5)]

        return {
            'p': d[0][:3], 'v': d[1][:3], 'a': d[2][:3],
            'j': d[3][:3], 's': d[4][:3],                   # jerk e snap: servono a DFL e controllo geometrico
            'psi': d[0][3], 'dpsi': d[1][3], 'ddpsi': d[2][3],
        }


class Hover:
    # T1: posizione e yaw costanti, stessa interfaccia di TraiettoriaInterpolata

    def __init__(self, p, psi:float = 0.0):
        self.p = np.array(p, dtype=float)
        self.psi = psi

        return

    def __call__(self, t:float):
        zero = np.zeros(3)
        return {
            'p': self.p.copy(), 'v': zero, 'a': zero.copy(),
            'j': zero.copy(), 's': zero.copy(),
            'psi': self.psi, 'dpsi': 0.0, 'ddpsi': 0.0,
        }


class Otto:
    # T2: lemniscata di Gerono nel piano orizzontale, con derivate analitiche fino allo snap
    #   x = cx + A sin(w t),   y = cy + (A/2) sin(2 w t),   z = z0,   psi costante
    # Il riferimento parte già in moto: a t = 0 è al centro con velocità (A w, A w, 0)

    def __init__(self, A:float, omega:float, centro = (0.0, 0.0, 1.5), psi:float = 0.0):
        self.A = A                                      # semiampiezza lungo x [m]
        self.omega = omega                              # [rad/s]
        self.centro = np.array(centro, dtype=float)     # (cx, cy, z0) [m]
        self.psi = psi
        self.periodo = 2*np.pi/omega                    # durata di un giro [s]

        return

    def __call__(self, t:float):
        A, w = self.A, self.omega
        s1, c1 = np.sin(w*t), np.cos(w*t)
        s2, c2 = np.sin(2*w*t), np.cos(2*w*t)

        # Derivate di ordine 0..4 di x e y
        x = [A*s1, A*w*c1, -A*w**2*s1, -A*w**3*c1, A*w**4*s1]
        y = [A/2*s2, A*w*c2, -2*A*w**2*s2, -4*A*w**3*c2, 8*A*w**4*s2]
        d = [np.array([x[r], y[r], 0.0]) for r in range(5)]

        return {
            'p': d[0] + self.centro, 'v': d[1], 'a': d[2],
            'j': d[3], 's': d[4],
            'psi': self.psi, 'dpsi': 0.0, 'ddpsi': 0.0,
        }


if __name__ == "__main__":
    # Confronto dei metodi sugli stessi punti: python comune/traiettoria.py
    import matplotlib.pyplot as plt

    punti = [[0, 0, 0], [1, 0, 1], [1, 1, 1.5], [0, 1, 1], [0, 0, 1]]
    tempi = [0, 2, 4, 6, 8]
    t = np.linspace(-0.5, 8.5, 2000)

    fig, ax = plt.subplots(4, 1, figsize=(9, 10), sharex=True)
    for metodo in METODI:
        traiettoria = TraiettoriaInterpolata(punti, metodo, tempi=tempi)
        rif = [traiettoria(ti) for ti in t]
        for i, chiave in enumerate(['p', 'v', 'a', 'j']):
            ax[i].plot(t, [r[chiave][0] for r in rif], label=metodo)

    for a, nome in zip(ax, ['x [m]', 'v_x [m/s]', 'a_x [m/s^2]', 'jerk_x [m/s^3]']):
        a.set_ylabel(nome)
        a.grid(True)
    ax[0].plot(tempi, [p[0] for p in punti], 'ko')
    ax[0].legend()
    ax[2].set_ylim(-3, 3)   # la lineare ha impulsi nei punti, qui si vedono come zeri
    ax[3].set_ylim(-5, 5)
    ax[-1].set_xlabel('t [s]')
    fig.tight_layout()
    plt.show()
