"""
Analisi lineare in hover (FASE 0 di STEPS_progetto.md), senza CoppeliaSim.

Modello: ogni canale è una catena di integratori con il motore del primo ordine. Gli ingressi sono
normalizzati in accelerazioni (u = tau/I oppure dT/m): i tre controllori cancellano massa e inerzia,
quindi qui non compaiono. L'analisi dipende solo da g, dai poli, da tau_m, dal passo T_c e da tau_f.

Convenzioni: z verso l'alto, angoli ZYX, hover con psi = 0, stati come scarti dall'hover.
Il controllo geometrico linearizzato coincide con il BS (stessa cascata, stesse derivate dell'assetto
desiderato), quindi non ha funzioni proprie: vale tutto quello che vale per il BS.

Uso:  python analisi/lineare.py   -> stampa i controlli contro i blocchi Verifica di STEPS_progetto.md
"""
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

# Le funzioni dai poli ai guadagni stanno in comune/taratura.py (FASE 1): qui si importano, non si riscrivono
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from comune.taratura import polinomio, zeta_min_bs, coefficienti_dfl, scala_poli
from comune.taratura import canale_bs as guadagni_bs

G = 9.81

# --- Parametri nominali dell'impianto (gli stessi di main_hover.py) ---
TAU_M = 0.02    # costante di tempo dei motori [s]
T_C = 0.005     # passo di controllo e di simulazione [s]
TAU_F = 0.02    # costante di tempo del filtro sulle derivate di phi_d, theta_d [s]

# --- Poli della FASE 0: (omega_n [rad/s], zeta) di ogni canale ---
# Sono quelli dei guadagni di main_hover.py: per esempio (k1, k2) = (1, 1.5) dà s^2 + 2.5 s + 2.5,
# cioè omega_n = sqrt(2.5) e zeta = 2.5/(2 sqrt(2.5)). Restano qui perché il rapporto della FASE 0 è scritto
# con la coppia del DFL su z a -9 doppio. La taratura in uso è comune/taratura.POLI (DFL su z a -20 doppio)
POLI = {
    'xy':    (np.sqrt(2.5), 2.5/(2*np.sqrt(2.5))),
    'z':     (np.sqrt(7.0), 5.0/(2*np.sqrt(7.0))),
    'att':   (9.0, 1.0),     # phi, theta
    'psi':   (5.0, 1.0),
    'dfl_z': (9.0, 1.0),     # coppia di poli in più del DFL su z (decisione aperta)
}

# Griglia di frequenze per margini e carta di Nichols [rad/s] (sotto la frequenza di Nyquist pi/T_C)
W = np.logspace(-2, np.log10(600), 4000)


# =====================================================================================
# Modello lineare dell'impianto
# =====================================================================================

def catena(tipo):
    """
    Parte rigida di un canale:  x' = A_r x + B_r a,  con a l'accelerazione realizzata dal motore.
    'laterale': stati (x, vx, theta, q), x'' = g theta, theta'' = a   (asse y uguale, con -g e phi)
    'doppio':   stati (q, dq), q'' = a                                (quota, imbardata, assetto da solo)
    """
    if tipo == 'laterale':
        A = np.zeros((4, 4))
        A[0, 1] = 1
        A[1, 2] = G
        A[2, 3] = 1
        B = np.array([0, 0, 0, 1.0])
    elif tipo == 'doppio':
        A = np.array([[0, 1.0], [0, 0]])
        B = np.array([0, 1.0])
    else:
        raise ValueError(tipo)
    return A, B


def impianto_continuo(tipo, tau_m=TAU_M):
    # Catena + motore a' = (u - a)/tau_m. Restituisce A, B e C_mis (stati misurati: tutti tranne a)
    A_r, B_r = catena(tipo)
    n = len(A_r)
    if tau_m == 0:
        return A_r, B_r, np.eye(n)

    A = np.zeros((n+1, n+1))
    A[:n, :n] = A_r
    A[:n, n] = B_r
    A[n, n] = -1/tau_m
    B = np.zeros(n+1)
    B[n] = 1/tau_m
    C_mis = np.hstack([np.eye(n), np.zeros((n, 1))])
    return A, B, C_mis


def impianto_discreto(tipo, T=T_C, tau_m=TAU_M):
    """
    Modello esatto a tempo discreto, nell'ordine del ciclo di main_hover.py e di Drone.comanda_motori:
      1. si legge lo stato a t_k  2. si calcola u_k  3. il motore si aggiorna: a_{k+1} = (1-al) a_k + al u_k
      4. la forza con il motore GIÀ aggiornato resta costante per tutto il passo.
    """
    A_r, B_r = catena(tipo)
    n = len(A_r)

    # Parte rigida con ingresso costante sul passo: esponenziale della matrice aumentata
    M = np.zeros((n+1, n+1))
    M[:n, :n] = A_r
    M[:n, n] = B_r
    E = expm(M*T)
    Phi, Gam = E[:n, :n], E[:n, n]

    if tau_m == 0:
        return Phi, Gam, np.eye(n)

    al = 1 - np.exp(-T/tau_m)       # come Drone.alpha
    F = np.zeros((n+1, n+1))
    F[:n, :n] = Phi
    F[:n, n] = (1 - al)*Gam
    F[n, n] = 1 - al
    H = np.zeros(n+1)
    H[:n] = al*Gam
    H[n] = al
    C_mis = np.hstack([np.eye(n), np.zeros((n, 1))])
    return F, H, C_mis


# =====================================================================================
# Leggi di controllo linearizzate
# =====================================================================================

class Legge:
    """
    Legge di controllo lineare dagli stati misurati y all'ingresso normalizzato u.
      continua:  xi' = Ac xi + Bc y        discreta:  xi[k+1] = Ac xi[k] + Bc y[k]
                 u   = Cc xi + Dc y                    u[k]    = Cc xi[k] + Dc y[k]
    Una legge statica ha solo Dc. 'disc' è la versione discreta, scritta come nel codice.
    """

    def __init__(self, Dc, Ac=None, Bc=None, Cc=None, xi0=None, disc=None):
        self.Dc = np.array(Dc, dtype=float)
        ny = len(self.Dc)
        if Ac is None:
            self.Ac, self.Bc, self.Cc = np.zeros((0, 0)), np.zeros((0, ny)), np.zeros(0)
        else:
            self.Ac, self.Bc, self.Cc = np.array(Ac, float), np.array(Bc, float), np.array(Cc, float)
        self.n = len(self.Ac)
        self._xi0 = xi0
        self.disc = disc if disc is not None else self    # le leggi statiche sono uguali nei due casi

    def xi0(self, y0):
        # Stato iniziale del controllore, dato lo stato misurato iniziale
        if self._xi0 is None:
            return np.zeros(self.n)
        return np.array(self._xi0(y0), float)

    def K(self, s):
        # Funzione di trasferimento da y a u su un vettore di s (continua) o di z (discreta): una riga per s
        s = np.atleast_1d(s)
        if self.n == 0:
            return np.tile(self.Dc.astype(complex), (len(s), 1))
        X = np.linalg.solve(s[:, None, None]*np.eye(self.n) - self.Ac, np.broadcast_to(self.Bc, (len(s), *self.Bc.shape)))
        return np.einsum('j,kjl->kl', self.Cc, X) + self.Dc


def legge_bs_laterale(poli=POLI, derivate='filtrate', tau_f=TAU_F, T=T_C):
    """
    BS sull'asse x linearizzato. y = (x, vx, theta, q).
      posizione:  mu = -alpha_p x - beta_p vx        inversione:  theta_d = mu/g
      assetto:    u = theta_d'' + alpha_a (theta_d - theta) + beta_a (theta_d' - q)
    derivate: 'nulle' (tau_filtro=None), 'esatte' (lungo il modello), 'filtrate' (come nel codice)
    """
    bp, ap = polinomio(*poli['xy'])
    ba, aa = polinomio(*poli['att'])
    Kd = np.array([-ap/G, -bp/G, 0, 0])     # theta_d = Kd y
    e_th = np.array([0, 0, 1.0, 0])
    e_q = np.array([0, 0, 0, 1.0])

    if derivate == 'nulle':
        return Legge(aa*(Kd - e_th) - ba*e_q)

    if derivate == 'esatte':
        dKd = np.array([0, -ap/G, -bp, 0])      # theta_d'  = -(alpha_p vx + beta_p g theta)/g
        ddKd = np.array([0, 0, -ap, -bp])       # theta_d'' = -(alpha_p theta + beta_p q)
        return Legge(ddKd + aa*(Kd - e_th) + ba*(dKd - e_q))

    if derivate != 'filtrate':
        raise ValueError(derivate)

    # Continua: due filtri s/(tau_f s + 1) in serie, stati xi = (w1, w2)
    #   d1 = (theta_d - w1)/tau_f,  w1' = d1        d2 = (d1 - w2)/tau_f,  w2' = d2
    tf = tau_f
    Ac = [[-1/tf, 0], [-1/tf**2, -1/tf]]
    Bc = [Kd/tf, Kd/tf**2]
    Cc = [-1/tf**2 - ba/tf, -1/tf]
    Dc = Kd/tf**2 + aa*(Kd - e_th) + ba*(Kd/tf - e_q)
    xi0_c = lambda y0: [Kd @ y0, 0.0]           # filtro già a regime su theta_d(0): nessun impulso iniziale

    # Discreta, come ControlloreBackstepping._derivate_assetto (Eulero all'indietro)
    #   d1 = a d1_prec + b (theta_d - theta_d_prec)     d2 = a d2_prec + b (d1 - d1_prec)
    # stati xi = (theta_d_prec, d1_prec, d2_prec); primo_passo: theta_d_prec = theta_d(0), derivate nulle
    a = tf/(tf + T)
    b = 1/(tf + T)
    Ad = [[0, 0, 0], [-b, a, 0], [-b**2, b*(a - 1), a]]
    Bd = [Kd, b*Kd, b**2*Kd]
    Cd = [-b**2 - ba*b, b*(a - 1) + ba*a, a]
    Dd = b**2*Kd + aa*(Kd - e_th) + ba*(b*Kd - e_q)
    disc = Legge(Dd, Ad, Bd, Cd, xi0=lambda y0: [Kd @ y0, 0.0, 0.0])

    return Legge(Dc, Ac, Bc, Cc, xi0=xi0_c, disc=disc)


def legge_dfl_laterale(poli=POLI):
    # DFL sull'asse x: x'''' = -c3 x''' - c2 x'' - c1 x' - c0 x, con x'' = g theta e x''' = g q dal modello
    c3, c2, c1, c0 = coefficienti_dfl(polinomio(*poli['xy']), polinomio(*poli['att']))
    return Legge([-c0/G, -c1/G, -c2, -c3])


def legge_doppio(pol):
    # BS (e GEO) su quota, imbardata o assetto da solo: u = -alpha q - beta dq. y = (q, dq)
    beta, alpha = polinomio(*pol)
    return Legge([-alpha, -beta])


def legge_dfl_z(poli=POLI, T=T_C):
    """
    DFL sulla quota, con il compensatore zeta = (T, T') normalizzato per m: xi = (dT/m, T'/m).
      xi1' = xi2,  xi2' = -c0 z - c1 vz - c2 xi1 - c3 xi2,  u = xi1 (spinta comandata)
    Il compensatore integra il comando e non vede il motore: z'' e z''' li ricava da xi.
    Discreta: Eulero in avanti, come in DFL_control/controllore_dfl.py. zeta(0) = (mg, 0) -> xi(0) = 0.
    """
    c3, c2, c1, c0 = coefficienti_dfl(polinomio(*poli['z']), polinomio(*poli['dfl_z']))
    Ac = np.array([[0, 1.0], [-c2, -c3]])
    Bc = np.array([[0, 0], [-c0, -c1]])
    Cc = np.array([1.0, 0])
    disc = Legge([0, 0], np.eye(2) + T*Ac, T*Bc, Cc)
    return Legge([0, 0], Ac, Bc, Cc, disc=disc)


# =====================================================================================
# Anello chiuso: poli e risposte
# =====================================================================================

def matrice_anello_chiuso(F, H, C_mis, legge):
    # Stato complessivo (impianto, controllore). Vale in continuo (A, B) e in discreto (Phi, Gamma)
    n, nc = len(F), legge.n
    M = np.zeros((n + nc, n + nc))
    M[:n, :n] = F + np.outer(H, legge.Dc @ C_mis)
    if nc:
        M[:n, n:] = np.outer(H, legge.Cc)
        M[n:, :n] = legge.Bc @ C_mis
        M[n:, n:] = legge.Ac
    return M


def poli_continui(tipo, legge, tau_m=TAU_M):
    A, B, C = impianto_continuo(tipo, tau_m)
    return ordina(np.linalg.eigvals(matrice_anello_chiuso(A, B, C, legge)))


def poli_discreti(tipo, legge, T=T_C, tau_m=TAU_M):
    # Autovalori z dell'anello discreto e loro equivalenti in s = ln(z)/T (gli z nulli si scartano)
    F, H, C = impianto_discreto(tipo, T, tau_m)
    z = np.linalg.eigvals(matrice_anello_chiuso(F, H, C, legge.disc))
    z_ok = z[np.abs(z) > 1e-9]
    return ordina(np.log(z_ok.astype(complex))/T), np.max(np.abs(z))


def ordina(p):
    # Dal polo più veloce al più lento; i complessi coniugati restano vicini
    p = np.array(p, complex)
    return p[np.lexsort((-p.imag, p.real))]


def risposta(tipo, legge, x0_rigido, durata, discreta=True, T=T_C, tau_m=TAU_M):
    """
    Risposta libera da uno scarto iniziale (equivale a un gradino del riferimento).
    discreta=True: anello esatto come nel ciclo in CoppeliaSim; False: anello continuo (senza campionamento).
    Restituisce t, stati dell'impianto (senza il motore) e comando u.
    """
    if discreta:
        F, H, C = impianto_discreto(tipo, T, tau_m)
        lg = legge.disc
        M = matrice_anello_chiuso(F, H, C, lg)
    else:
        A, B, C = impianto_continuo(tipo, tau_m)
        lg = legge
        M = expm(matrice_anello_chiuso(A, B, C, lg)*T)

    n = len(C[0])
    x0 = np.zeros(n)
    x0[:len(x0_rigido)] = x0_rigido             # motore all'hover: a(0) = 0
    stato = np.concatenate([x0, lg.xi0(C @ x0)])

    N = int(round(durata/T)) + 1
    X = np.zeros((N, len(stato)))
    X[0] = stato
    for k in range(N - 1):
        X[k+1] = M @ X[k]

    Y = X[:, :n] @ C.T
    U = Y @ lg.Dc + (X[:, n:] @ lg.Cc if lg.n else 0)
    return np.arange(N)*T, Y, U


def tempo_entro(t, e, soglia):
    # Primo istante dopo il quale |e| resta sotto la soglia
    fuori = np.nonzero(np.abs(e) > soglia)[0]
    if len(fuori) == 0:
        return 0.0
    if fuori[-1] == len(e) - 1:
        return np.inf
    return t[fuori[-1] + 1]


# =====================================================================================
# Risposta in frequenza e margini
# =====================================================================================

def _impianto_in_frequenza(F, H, C, s):
    # y(s) = C (sI - F)^-1 H per ogni s: una riga per s
    X = np.linalg.solve(s[:, None, None]*np.eye(len(F)) - F, np.broadcast_to(H[:, None], (len(s), len(H), 1)))
    return X[:, :, 0] @ C.T


def anello_attuatore(tipo, legge, w=W, tau_m=TAU_M, T=T_C, ritardo=True):
    # Anello aperto all'ingresso dell'impianto: L(jw) = -K(jw) C (jwI - A)^-1 B e^{-jw T/2}
    A, B, C = impianto_continuo(tipo, tau_m)
    s = 1j*np.asarray(w)
    L = -np.sum(legge.K(s)*_impianto_in_frequenza(A, B, C, s), axis=1)
    if ritardo:
        L = L*np.exp(-s*T/2)
    return L


def anello_attuatore_discreto(tipo, legge, w=None, T=T_C, tau_m=TAU_M):
    # Lo stesso taglio sul modello discreto esatto, con z = e^{jwT}, fino alla frequenza di Nyquist pi/T
    if w is None:
        w = np.logspace(-2, np.log10(np.pi/T*(1 - 1e-9)), 4000)
    F, H, C = impianto_discreto(tipo, T, tau_m)
    z = np.exp(1j*np.asarray(w)*T)
    return -np.sum(legge.disc.K(z)*_impianto_in_frequenza(F, H, C, z), axis=1)


def anello_posizione_bs(poli=POLI, derivate='filtrate', w=W, tau_f=TAU_F, tau_m=TAU_M, T=T_C):
    """
    Anello di posizione del BS aperto dove entra theta_d:  L = (beta_p s + alpha_p)/s^2 * T_att(s)
    T_att = N G / (1 + (beta_a s + alpha_a) G),  G = e^{-sT/2}/(s^2 (tau_m s + 1))
    N = alpha_a con derivate nulle, N = D^2 + beta_a D + alpha_a con D = s/(tau_f s + 1) filtrate
    """
    bp, ap = polinomio(*poli['xy'])
    ba, aa = polinomio(*poli['att'])
    s = 1j*w
    Gs = np.exp(-s*T/2)/(s**2*(tau_m*s + 1))
    if derivate == 'nulle':
        N = aa
    elif derivate == 'filtrate':
        D = s/(tau_f*s + 1)
        N = D**2 + ba*D + aa
    else:
        raise ValueError(derivate)
    T_att = N*Gs/(1 + (ba*s + aa)*Gs)
    return (bp*s + ap)/s**2*T_att


def margini(L, w=W):
    """
    Margine di fase [deg] alla pulsazione di attraversamento (il minimo se ce n'è più d'una),
    margine di guadagno superiore [dB] e, per gli anelli condizionatamente stabili, quello inferiore
    (di quanto si può ridurre il guadagno, in dB negativi).
    """
    mod = 20*np.log10(np.abs(L))
    fase = np.degrees(np.unwrap(np.angle(L)))

    pm, wc = np.inf, None
    for i in np.nonzero(np.diff(np.sign(mod)) != 0)[0]:
        f = mod[i]/(mod[i] - mod[i+1])
        w_i = w[i] + f*(w[i+1] - w[i])
        fase_i = fase[i] + f*(fase[i+1] - fase[i])
        pm_i = (fase_i + 180 + 180) % 360 - 180
        if pm_i < pm:
            pm, wc = pm_i, w_i

    gm_sup, wg, gm_inf = np.inf, None, None
    giro = np.floor((fase + 180)/360)
    for i in np.nonzero(np.diff(giro) != 0)[0]:
        f_cr = 360*max(giro[i], giro[i+1]) - 180
        f = (f_cr - fase[i])/(fase[i+1] - fase[i])
        g_i = -(mod[i] + f*(mod[i+1] - mod[i]))
        if g_i > 0 and g_i < gm_sup:
            gm_sup, wg = g_i, w[i] + f*(w[i+1] - w[i])
        if g_i < 0 and (gm_inf is None or g_i > gm_inf):
            gm_inf = g_i

    return {'pm': pm, 'wc': wc, 'gm': gm_sup, 'wg': wg, 'gm_inf': gm_inf}


def curva_nichols(L, w=W):
    # Fase in gradi, spostata di multipli di 360 perché alla pulsazione di attraversamento stia in (-360, 0]
    mod = 20*np.log10(np.abs(L))
    fase = np.degrees(np.unwrap(np.angle(L)))
    i_c = np.argmin(np.abs(mod))
    fase = fase - 360*np.ceil(fase[i_c]/360)
    return fase, mod


# =====================================================================================
# Le analisi della FASE 0
# =====================================================================================

def pm_assetto(wn, zeta=1.0, tau_m=TAU_M, T=T_C):
    return margini(anello_attuatore('doppio', legge_doppio((wn, zeta)), tau_m=tau_m, T=T))


def wn_max_assetto(pm_min=50.0, zeta=1.0):
    # Esercizio 8: la omega_n più alta dell'assetto con margine di fase >= pm_min (bisezione)
    a, b = 3.0, 30.0
    for _ in range(40):
        c = 0.5*(a + b)
        if pm_assetto(c, zeta)['pm'] >= pm_min:
            a = c
        else:
            b = c
    return 0.5*(a + b)


def tau_m_instabile(poli=POLI):
    # Esercizio 13: tau_m oltre il quale l'asse x (BS, derivata filtrata) diventa instabile
    def raggio(tm):
        return poli_discreti('laterale', legge_bs_laterale(poli), tau_m=tm)[1]
    def pm_act(tm):
        return margini(anello_attuatore('laterale', legge_bs_laterale(poli), tau_m=tm))['pm']
    a, b = 0.02, 0.4
    for _ in range(40):
        c = 0.5*(a + b)
        a, b = (c, b) if raggio(c) < 1 else (a, c)
    tm_disc = 0.5*(a + b)
    a, b = 0.02, 0.4
    for _ in range(40):
        c = 0.5*(a + b)
        a, b = (c, b) if pm_act(c) > 0 else (a, c)
    return tm_disc, 0.5*(a + b)


def tabella_guadagni(poli=POLI):
    righe = {}
    for canale in ['xy', 'z', 'att', 'psi']:
        wn, zeta = poli[canale]
        beta, alpha = polinomio(wn, zeta)
        righe[canale] = {'wn': wn, 'zeta': zeta, 'beta': beta, 'alpha': alpha,
                         'zeta_min': zeta_min_bs(wn), 'k_bs': guadagni_bs(wn, zeta),
                         'poli': ordina(np.roots([1, beta, alpha]))}
    righe['dfl_xy'] = coefficienti_dfl(polinomio(*poli['xy']), polinomio(*poli['att']))
    righe['dfl_z'] = coefficienti_dfl(polinomio(*poli['z']), polinomio(*poli['dfl_z']))
    righe['dfl_psi'] = np.array(polinomio(*poli['psi']))
    return righe


def gradino_x(legge, discreta=True, tau_m=TAU_M, ampiezza=0.1, durata=10.0):
    t, Y, U = risposta('laterale', legge, [ampiezza, 0, 0, 0], durata, discreta, tau_m=tau_m)
    return t, Y[:, 0], Y[:, 2], U


def gradino_z(legge, discreta=True, ampiezza=0.5, durata=4.0):
    t, Y, U = risposta('doppio', legge, [ampiezza, 0], durata, discreta)
    return t, Y[:, 0], U


def confronto_derivate(poli=POLI):
    # Esercizio 11: BS sull'asse x con derivate esatte, filtrate, nulle (continuo, con il motore)
    esito = {}
    t, x_es, _, _ = gradino_x(legge_bs_laterale(poli, 'esatte'), discreta=False)
    for d in ['esatte', 'filtrate', 'nulle']:
        lg = legge_bs_laterale(poli, d)
        t, x, th, _ = gradino_x(lg, discreta=False)
        esito[d] = {'poli_senza_motore': poli_continui('laterale', lg, tau_m=0),
                    'poli_con_motore': poli_continui('laterale', lg),
                    'dist_max': np.max(np.abs(x - x_es)),
                    't_5mm': tempo_entro(t, x, 0.005)}
    return esito


def confronto_dfl_z(poli=POLI, extra=(9.0, 15.0, 25.0)):
    # Esercizio 12: gradino di 0.5 m su z, BS/GEO contro DFL con diverse coppie di poli in più
    t, z_bs, _ = gradino_z(legge_doppio(poli['z']), discreta=False)
    esito = {'bs': {'t_5': tempo_entro(t, z_bs, 0.025)}}
    for lam in extra:
        p = dict(poli)
        p['dfl_z'] = (lam, 1.0)
        t, z, _ = gradino_z(legge_dfl_z(p), discreta=False)
        esito[lam] = {'t_5': tempo_entro(t, z, 0.025), 'dist_max': np.max(np.abs(z - z_bs))}
    return esito


def margini_tau_m(poli=POLI, valori=(0.02, 0.04, 0.06, 0.08)):
    # Esercizio 13: margini al variare della costante di tempo dei motori
    esito = {}
    for tm in valori:
        esito[tm] = {'att': margini(anello_attuatore('doppio', legge_doppio(poli['att']), tau_m=tm)),
                     'act': margini(anello_attuatore('laterale', legge_bs_laterale(poli), tau_m=tm)),
                     'act_dfl': margini(anello_attuatore('laterale', legge_dfl_laterale(poli), tau_m=tm)),
                     'z': margini(anello_attuatore('doppio', legge_doppio(poli['z']), tau_m=tm)),
                     'z_dfl': margini(anello_attuatore('doppio', legge_dfl_z(poli), tau_m=tm)),
                     'psi': margini(anello_attuatore('doppio', legge_doppio(poli['psi']), tau_m=tm))}
    return esito


def sensibilita_kappa(valori=(0.7, 1.0, 1.3)):
    # Esercizio 14: poli scalati di kappa
    esito = {}
    for k in valori:
        p = scala_poli(POLI, k)
        lg = legge_bs_laterale(p)
        t, x, th, _ = gradino_x(lg)
        esito[k] = {'poli': p, 'guadagni': tabella_guadagni(p),
                    'pm_att': margini(anello_attuatore('doppio', legge_doppio(p['att'])))['pm'],
                    'pm_act': margini(anello_attuatore('laterale', lg))['pm'],
                    'pm_pos': margini(anello_posizione_bs(p))['pm'],
                    't_5mm': tempo_entro(t, x, 0.005),
                    'incl_max': np.degrees(np.max(np.abs(th)))}
    return esito


def massa_non_nota(fattori=(1.5, 2.0, 2.5), poli=POLI):
    """
    Controprova (planning §9.4): il controllore usa massa e inerzia nominali, il drone pesa k volte tanto.
    - quota: nessuno ha integrale, resta un errore a regime (il drone sta sotto il riferimento)
        BS/GEO:  m0 (alpha_z e + g) = m g          ->  e = g (k - 1)/alpha_z
        DFL:     c2 (z''_stimata) + c0 e = 0, con z''_stimata = g (k - 1)  ->  e = (c2/c0) g (k - 1)
    - assetto: tau = J0 nu, quindi il guadagno d'anello si riduce di 1/k
    """
    _, az = polinomio(*poli['z'])
    c3, c2, c1, c0 = coefficienti_dfl(polinomio(*poli['z']), polinomio(*poli['dfl_z']))
    L_att = anello_attuatore('doppio', legge_doppio(poli['att']))
    L_act = anello_attuatore('laterale', legge_bs_laterale(poli))
    esito = {}
    for k in fattori:
        esito[k] = {'ez_bs': G*(k - 1)/az, 'ez_dfl': (c2/c0)*G*(k - 1),
                    'att': margini(L_att/k), 'act': margini(L_act/k)}
    return esito


def controlli():
    """
    I numeri dei blocchi Verifica di STEPS_progetto.md, confrontati con quelli calcolati.
    Ogni riga: (esercizio, descrizione, atteso, calcolato, tolleranza assoluta)
    """
    r = []
    gu = tabella_guadagni()
    r.append(('7', 'zeta minima del BS a omega_n = 9', 0.994, zeta_min_bs(9.0), 0.001))
    r.append(('7', 'zeta minima del BS a omega_n = 1.58', 0.775, zeta_min_bs(1.58), 0.001))

    m9 = pm_assetto(9.0)
    r.append(('8', 'omega_n massima dell\'assetto con PM >= 50° [rad/s]', 10.8, wn_max_assetto(), 0.05))
    r.append(('8', 'PM dell\'assetto a omega_n = 9 [deg]', 54, m9['pm'], 0.5))
    r.append(('8', 'GM dell\'assetto a omega_n = 9 [dB]', 26, m9['gm'], 0.5))

    r.append(('9', 'PM di posizione su theta_d, derivate nulle [deg]', 37, margini(anello_posizione_bs(derivate='nulle'))['pm'], 0.5))
    r.append(('9', 'PM di posizione su theta_d, filtrate 20 ms [deg]', 70, margini(anello_posizione_bs())['pm'], 1.0))
    r.append(('9', 'PM all\'attuatore, filtrate 20 ms [deg]', 43, margini(anello_attuatore('laterale', legge_bs_laterale()))['pm'], 0.5))
    r.append(('9', 'PM all\'attuatore, filtrate 10 ms [deg] (HTML)', 45, margini(anello_attuatore('laterale', legge_bs_laterale(tau_f=0.01)))['pm'], 0.5))

    for nome, att, calc in zip(['k1 x', 'k2 x', 'k1 z', 'k2 z', 'k1 att', 'k2 att', 'k1 psi', 'k2 psi'],
                               [1, 1.5, 2, 3, 8, 10, 4, 6],
                               [*gu['xy']['k_bs'], *gu['z']['k_bs'], *gu['att']['k_bs'], *gu['psi']['k_bs']]):
        r.append(('10', 'BS '+nome, att, calc, 1e-6))
    for i, att in enumerate([20.5, 128.5, 247.5, 202.5]):
        r.append(('10', 'DFL x: c'+str(3-i), att, gu['dfl_xy'][i], 1e-6))

    der = confronto_derivate()
    pn = der['nulle']['poli_senza_motore']
    pr = np.sort(pn[np.abs(pn.imag) < 1e-6].real)
    r.append(('11', 'derivate nulle: polo reale veloce', -12.8, pr[0], 0.05))
    r.append(('11', 'derivate nulle: polo reale lento', -1.9, pr[-1], 0.05))
    pc = pn[np.abs(pn.imag) > 1e-6]
    r.append(('11', 'derivate nulle: parte reale dei complessi', -1.64, pc[0].real, 0.01))
    r.append(('11', 'derivate nulle: parte immaginaria dei complessi', 2.37, abs(pc[0].imag), 0.01))
    r.append(('11', 'distanza massima, filtrate [mm]', 0.4, 1000*der['filtrate']['dist_max'], 0.1))
    r.append(('11', 'distanza massima, nulle [mm]', 19, 1000*der['nulle']['dist_max'], 0.5))

    dz = confronto_dfl_z()
    r.append(('12', 'BS/GEO z: entro il 5% [s]', 1.64, dz['bs']['t_5'], 0.01))
    for lam, t_at, d_at in [(9.0, 1.88, 109), (15.0, 1.78, 70), (25.0, 1.72, 43)]:
        r.append(('12', 'DFL z poli -'+str(int(lam))+': entro il 5% [s]', t_at, dz[lam]['t_5'], 0.01))
        r.append(('12', 'DFL z poli -'+str(int(lam))+': distanza da BS [mm]', d_at, 1000*dz[lam]['dist_max'], 1))

    mt = margini_tau_m()
    for tm, act, att in [(0.02, 43, 54), (0.04, 28, 40), (0.06, 18, 30), (0.08, 11, 23)]:
        r.append(('13', 'tau_m = '+str(int(tm*1000))+' ms: PM all\'attuatore [deg]', act, mt[tm]['act']['pm'], 0.5))
        r.append(('13', 'tau_m = '+str(int(tm*1000))+' ms: PM dell\'assetto [deg]', att, mt[tm]['att']['pm'], 0.5))
    r.append(('13', 'tau_m di instabilità (PM = 0, continuo) [ms]', 120, 1000*tau_m_instabile()[1], 10))

    sk = sensibilita_kappa()
    for k, att, act, pos, t5, inc in [(0.7, 60, 51, 70, 3.38, 0.4), (1.0, 54, 43, 70, 2.37, 0.8), (1.3, 48, 36, 73, 2.07, 1.4)]:
        r.append(('14', 'kappa '+str(k)+': PM assetto [deg]', att, sk[k]['pm_att'], 0.5))
        r.append(('14', 'kappa '+str(k)+': PM attuatore [deg]', act, sk[k]['pm_act'], 0.5))
        r.append(('14', 'kappa '+str(k)+': PM posizione [deg]', pos, sk[k]['pm_pos'], 1.0))
        r.append(('14', 'kappa '+str(k)+': entro 5 mm [s]', t5, sk[k]['t_5mm'], 0.02))
        r.append(('14', 'kappa '+str(k)+': inclinazione massima [deg]', inc, sk[k]['incl_max'], 0.06))   # 1.35 arrotondato a 1.4 nel Verifica

    return r


if __name__ == "__main__":
    np.set_printoptions(precision=3, suppress=True)

    print("\n--- Controlli contro i blocchi Verifica di STEPS_progetto.md ---")
    ok_tutti = True
    for es, descr, atteso, calc, tol in controlli():
        ok = abs(calc - atteso) <= tol
        ok_tutti &= ok
        print(("ok " if ok else "NO ")+"Es."+es.ljust(3)+descr.ljust(58)+" atteso "+str(atteso).rjust(7)+"  calcolato "+f"{calc:9.3f}")
    print("\nTutti i controlli passano" if ok_tutti else "\nALCUNI CONTROLLI NON PASSANO")

    print("\n--- Poli dell'asse x (BS, derivata filtrata) ---")
    lg = legge_bs_laterale()
    print("continuo, con motore:  "+str(poli_continui('laterale', lg)))
    s_d, raggio = poli_discreti('laterale', lg)
    print("discreto esatto (s = ln z / T):  "+str(s_d)+"   raggio spettrale "+f"{raggio:.4f}")

    print("\n--- Margini dell'asse x all'attuatore: continuo con mezzo passo di ritardo contro discreto esatto ---")
    print("continuo: "+str({k: round(v, 2) for k, v in margini(anello_attuatore('laterale', lg)).items() if v is not None}))
    print("discreto: "+str({k: round(v, 2) for k, v in margini(anello_attuatore_discreto('laterale', lg)).items() if v is not None}))
