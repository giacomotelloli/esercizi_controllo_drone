"""
Verifica offline delle tre classi (FASE 4.1 di STEPS_progetto.md), senza CoppeliaSim.

Le vere classi dei controllori, con guadagni e limiti presi da comune/taratura.py, comandano un impianto
non lineare in Python con la stessa catena del Drone: mixer a X, saturazione di Omega, motore del primo
ordine aggiornato prima di applicare la forza, corpo rigido integrato con RK4. Le risposte si confrontano
con le previsioni di analisi/lineare.py calcolate con gli stessi poli.

Massa e inerzia sono di prova: per il confronto con la previsione lineare non contano, perché i tre
controllori le cancellano.

Uso:  python analisi/verifica_controllori.py
"""
import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RADICE))
sys.path.insert(0, str(RADICE / 'analisi'))
sys.path.insert(0, str(RADICE / 'Backstepping_control'))
sys.path.insert(0, str(RADICE / 'DFL_control'))
sys.path.insert(0, str(RADICE / 'Geometric_control'))

import lineare as L
from comune import taratura
from comune.rotazioni import hat, matrice_R, angoli_da_R, angolo_in_pi
from comune.traiettoria import TraiettoriaInterpolata
from controllore_backstepping import ControlloreBackstepping
from controllore_dfl import ControlloreDFL
from controllore_geometrico import ControlloreGeometrico

G = taratura.G
DT = 0.005
CLASSI = {'BS': ControlloreBackstepping, 'DFL': ControlloreDFL, 'GEO': ControlloreGeometrico}

# Drone di prova
M_PROVA = 1.0
J_PROVA = (0.0082, 0.0085, 0.0149)


class Impianto:
    def __init__(self, m=M_PROVA, J=J_PROVA, a=0.12, k_f=1e-5, k_m=1.6e-7, Omega_max=1000.0, tau_m=0.02):
        self.m = m
        self.J = np.diag(J)
        self.J_inv = np.linalg.inv(self.J)
        pos = np.array([[a, a], [-a, a], [-a, -a], [a, -a]])
        segni = np.array([1, -1, 1, -1])
        self.M = np.array([k_f*np.ones(4), k_f*pos[:, 1], -k_f*pos[:, 0], k_m*segni])
        self.M_inv = np.linalg.inv(self.M)
        self.Omega_max = Omega_max
        self.T_max = 4*k_f*Omega_max**2
        self.alpha = 1 - np.exp(-DT/tau_m)
        self.reset(np.zeros(3), np.zeros(3))

    def reset(self, p0, eta0):
        self.p = np.array(p0, float)
        self.v = np.zeros(3)
        self.R = matrice_R(eta0)
        self.w = np.zeros(3)
        self.Omega = np.sqrt(np.maximum(self.M_inv @ np.array([self.m*G, 0, 0, 0]), 0))   # motori all'hover

    def stato(self):
        return {'p': self.p.copy(), 'v': self.v.copy(), 'R': self.R.copy(), 'eta': angoli_da_R(self.R), 'omega': self.w.copy()}

    def passo(self, T, tau):
        # Restituisce True se la richiesta del mixer esce da [0, Omega_max]
        O2 = np.maximum(self.M_inv @ np.array([T, *tau]), 0)
        Ocmd = np.clip(np.sqrt(O2), 0, self.Omega_max)
        self.Omega = self.Omega + self.alpha*(Ocmd - self.Omega)
        u = self.M @ self.Omega**2
        f, tq = u[0], u[1:]
        # corpo rigido con forza e coppia costanti sul passo: RK4 su 5 sottopassi
        h = DT/5
        for _ in range(5):
            def der(p, v, R, w):
                return v, -G*np.array([0, 0, 1.0]) + f/self.m*R[:, 2], R @ hat(w), self.J_inv @ (tq - np.cross(w, self.J @ w))
            k1 = der(self.p, self.v, self.R, self.w)
            k2 = der(*[x + h/2*k for x, k in zip((self.p, self.v, self.R, self.w), k1)])
            k3 = der(*[x + h/2*k for x, k in zip((self.p, self.v, self.R, self.w), k2)])
            k4 = der(*[x + h*k for x, k in zip((self.p, self.v, self.R, self.w), k3)])
            self.p, self.v, self.R, self.w = [x + h/6*(a + 2*b + 2*c + d) for x, a, b, c, d in zip((self.p, self.v, self.R, self.w), k1, k2, k3, k4)]
            U, _, Vt = np.linalg.svd(self.R)
            self.R = U @ Vt
        return np.max(Ocmd) >= self.Omega_max - 1e-9 or np.min(O2) <= 0


def crea(nome, imp, poli=taratura.POLI):
    # Il controllore come lo creerà esecuzione.py: tutto da taratura.parametri
    return CLASSI[nome](**taratura.parametri(nome, imp.m, imp.J, imp.T_max, DT, poli=poli))


def rif_fermo(p, psi=0.0):
    z = np.zeros(3)
    return lambda t: {'p': np.array(p, float), 'v': z, 'a': z, 'j': z, 's': z, 'psi': psi, 'dpsi': 0.0, 'ddpsi': 0.0}


def prova(nome, p0, eta0, rif, durata, imp=None, poli=taratura.POLI):
    imp = imp or Impianto()
    imp.reset(p0, eta0)
    c = crea(nome, imp, poli)
    N = int(round(durata/DT)) + 1
    P, E, S = np.zeros((N, 3)), np.zeros((N, 3)), np.zeros(N, bool)
    for k in range(N):
        t = k*DT
        s = imp.stato()
        P[k], E[k] = s['p'], s['eta']
        T, tau = c.calcola(s, rif(t))
        S[k] = imp.passo(T, tau)
        if not np.all(np.isfinite(imp.p)) or np.linalg.norm(imp.p) > 1e3:
            P[k:], E[k:] = np.nan, np.nan
            break
    return np.arange(N)*DT, P, E, S


def previsioni(poli=taratura.POLI, tau_f=taratura.TAU_F):
    # Risposte previste dal modello lineare discreto, con gli stessi poli delle classi
    bs_x = L.legge_bs_laterale(poli, tau_f=tau_f)
    t_x, x_bs, th_bs, _ = L.gradino_x(bs_x)
    _, x_dfl, th_dfl, _ = L.gradino_x(L.legge_dfl_laterale(poli))
    t_z, z_bs, _ = L.gradino_z(L.legge_doppio(poli['z']), durata=6.0)
    _, z_dfl, _ = L.gradino_z(L.legge_dfl_z(poli), durata=6.0)
    t_psi, Y_psi, _ = L.risposta('doppio', L.legge_doppio(poli['psi']), [np.radians(10), 0], 3.0)
    return {'x': {'t': t_x, 'BS': x_bs, 'DFL': x_dfl, 'GEO': x_bs, 'th_BS': th_bs, 'th_DFL': th_dfl},
            'z': {'t': t_z, 'BS': z_bs, 'DFL': z_dfl, 'GEO': z_bs},
            'psi': {'t': t_psi, 'BS': Y_psi[:, 0], 'DFL': Y_psi[:, 0], 'GEO': Y_psi[:, 0]}}


def gradini(prev=None):
    # Scarto massimo tra le classi e la propria previsione lineare: x, y (10 cm), z (0.5 m), yaw (10°)
    prev = prev or previsioni()
    rif0 = rif_fermo([0, 0, 0])
    esito = {}
    for nome in CLASSI:
        _, Px, Ex, _ = prova(nome, [0.1, 0, 0], [0, 0, 0], rif0, 10.0)
        _, Py, _, _ = prova(nome, [0, 0.1, 0], [0, 0, 0], rif0, 10.0)
        _, Pz, _, Sz = prova(nome, [0, 0, 0.5], [0, 0, 0], rif0, 6.0)
        _, _, Epsi, _ = prova(nome, [0, 0, 0], [0, 0, np.radians(10)], rif0, 3.0)
        esito[nome] = {'x_mm': 1000*np.max(np.abs(Px[:, 0] - prev['x'][nome])),
                       'y_mm': 1000*np.max(np.abs(Py[:, 1] - prev['x'][nome])),
                       'z_mm': 1000*np.max(np.abs(Pz[:, 2] - prev['z'][nome])),
                       'psi_deg': np.degrees(np.max(np.abs(Epsi[:, 2] - prev['psi'][nome]))),
                       'x_sim': Px[:, 0], 'z_sim': Pz[:, 2], 'sat_z': int(Sz.sum()),
                       'theta_max_deg': np.degrees(np.max(np.abs(Ex[:, 1])))}
    return esito


def t1_n2():
    # T1-N2: scostato di (1, 1, 0.5) m e 90° di yaw, riferimento fermo
    esito = {}
    for nome in CLASSI:
        t, P, E, S = prova(nome, [1, 1, 0.5], [0, 0, np.radians(90)], rif_fermo([0, 0, 0]), 10.0)
        esito[nome] = {'errore_finale_mm': 1000*np.linalg.norm(P[-1]), 'yaw_finale_deg': np.degrees(E[-1, 2]),
                       'incl_max_deg': np.degrees(np.max(np.abs(E[:, :2]))), 'passi_saturi': int(S.sum())}
    return esito


if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)
    rif0 = rif_fermo([0, 0, 0])

    print("\n=== 1-4. Gradini contro la previsione lineare discreta (stessi poli di comune/taratura.py) ===")
    for nome, e in gradini().items():
        print(f"{nome:4s} x {e['x_mm']:.4f} mm | y {e['y_mm']:.4f} mm | z {e['z_mm']:.3f} mm | yaw {e['psi_deg']:.4f} deg"
              f" | theta max {e['theta_max_deg']:.3f} deg | saturazione in z: {e['sat_z']} passi")

    print("\n=== 5. T1-N2: scostato di (1, 1, 0.5) m e 90 gradi di yaw ===")
    for nome, e in t1_n2().items():
        print(f"{nome:4s} errore a 10 s {e['errore_finale_mm']:.2f} mm, yaw {e['yaw_finale_deg']:.3f} deg"
              f" | incl max {e['incl_max_deg']:.1f} deg | passi in saturazione {e['passi_saturi']}")

    print("\n=== 6. Inseguimento: waypoint min_snap a 1.5 m/s (servono jerk e snap) ===")
    tr = TraiettoriaInterpolata([[0, 0, 0, 0], [2, 0, 1, 0.5], [2, 2, 1.5, 1.0], [0, 2, 1, 0.5], [0, 0, 0, 0]], 'min_snap', v_media=1.5, t_inizio=0.5)
    for nome in CLASSI:
        t, P, E, S = prova(nome, [0, 0, 0], [0, 0, 0], tr, tr.t_fine + 2)
        Pd = np.array([tr(ti)['p'] for ti in t]); psid = np.array([tr(ti)['psi'] for ti in t])
        err = np.linalg.norm(P - Pd, axis=1)
        print(f"{nome:4s} errore RMS {1000*np.sqrt(np.mean(err**2)):.2f} mm, max {1000*np.max(err):.2f} mm | yaw max {np.degrees(np.max(np.abs(angolo_in_pi(E[:,2]-psid)))):.3f} deg | incl max {np.degrees(np.max(np.abs(E[:,:2]))):.1f} deg")

    print("\n=== 7. T1-N3: rilasciato inclinato (rollio / beccheggio), riferimento fermo ===")
    for ang in [30, 60, 90, 120, 179]:
        riga = f"{ang:3d} deg:"
        for nome in CLASSI:
            for asse, eta0 in [('R', [np.radians(ang), 0, 0]), ('P', [0, np.radians(89.9 if ang == 90 else ang), 0])]:
                with np.errstate(all='ignore'):
                    try:
                        t, P, E, S = prova(nome, [0, 0, 0], eta0, rif0, 8.0)
                        ok = np.all(np.isfinite(P)) and np.linalg.norm(P[-1]) < 0.05
                        quota = np.nanmin(P[:, 2])
                        riga += f"  {nome}-{asse} {'ok ' if ok else 'NO '}(z min {quota:6.2f} m)"
                    except np.linalg.LinAlgError:
                        riga += f"  {nome}-{asse} singolare"
        print(riga)
