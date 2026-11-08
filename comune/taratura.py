"""
Taratura dei tre controllori (FASE 1 di STEPS_progetto.md): dai poli ai guadagni.

Il progetto si fa sul modello semplificato della slide 2 del PDF 09 (Controllers comparison): aerodinamica
ed effetti giroscopici trascurati, angoli piccoli (eta_dot ~ omega), drone simmetrico. In hover, con z verso l'alto:

    x'' = g theta      y'' = -g phi      z'' = dT/m
    phi'' = tau_phi/Ix      theta'' = tau_theta/Iy      psi'' = tau_psi/Iz

Ogni canale si chiude con lo stesso polinomio s^2 + beta s + alpha (beta = 2 zeta wn, alpha = wn^2) per BS,
DFL e GEO: in hover i tre diventano la stessa legge lineare, e le differenze nelle prove vengono solo dalla
costruzione dei metodi (STEPS 0.3). Unica eccezione la quota del DFL, che ha due poli in più.

Una sola fonte di verità: POLI. I guadagni si calcolano, non si scrivono mai a mano.

Uso, dopo aver creato il Drone (GEO ha bisogno di m e J della scena):
    from comune import taratura
    ctrl = ControlloreBackstepping(**taratura.parametri('BS', drone.massa, drone.J, drone.T_max, drone.dt))

    python -m comune.taratura   -> stampa poli e guadagni (GEO per unità di massa e d'inerzia)
"""
import numpy as np

G = 9.81

# --- Poli in anello chiuso: (omega_n [rad/s], zeta) di ogni canale, scelti con l'analisi lineare (FASE 0) ---
POLI = {
    'xy':    (np.sqrt(2.5), 2.5/(2*np.sqrt(2.5))),  # s^2 + 2.5 s + 2.5 -> -1.25 +- 0.97j (1.58 rad/s, zeta 0.79)
    'z':     (np.sqrt(7.0), 5.0/(2*np.sqrt(7.0))),  # s^2 + 5 s + 7     -> -2.5 +- 0.87j  (2.65 rad/s, zeta 0.945)
    'att':   (9.0, 1.0),                            # phi, theta: -9 doppio
    'psi':   (5.0, 1.0),                            # psi: -5 doppio
    'dfl_z': (20.0, 1.0),                           # coppia in più del DFL su z: -20 doppio, la più lenta con PM >= 50°
}

# --- Derivate dell'assetto desiderato (BS e GEO): derivata filtrata s/(tau_f s + 1) ---
TAU_F = 0.02    # [s], 4 passi di controllo

# --- Limiti, uguali per tutti ---
ANGOLO_MAX = np.radians(80)     # solo BS: inclinazione massima della spinta desiderata [rad]
FRAZIONE_T_MIN = 0.1            # spinta minima = 0.1 mg (nel BS è già il limite su F_z)


def polinomio(wn, zeta):
    # Coefficienti (beta, alpha) di s^2 + beta s + alpha
    return 2*zeta*wn, wn**2


def polinomi(poli=POLI):
    # (beta, alpha) di ogni canale
    return {canale: polinomio(*wn_zeta) for canale, wn_zeta in poli.items()}


def zeta_min_bs(wn):
    # Il canale del BS dà s^2 + (k1 + k2) s + (1 + k1 k2): k1, k2 reali e positivi solo se zeta >= sqrt(1 - 1/wn^2)
    if wn <= 1:
        return np.inf
    return np.sqrt(1 - 1/wn**2)


def canale_bs(wn, zeta):
    # (k1, k2) del BS su un canale: k1 + k2 = beta e 1 + k1 k2 = alpha, quindi radici di k^2 - beta k + (alpha - 1) = 0
    beta, alpha = polinomio(wn, zeta)
    delta = beta**2 - 4*(alpha - 1)
    if alpha <= 1 or delta < -1e-9:
        raise ValueError("Poli fuori dal vincolo del BS: wn="+str(round(wn,3))+" zeta="+str(round(zeta,3))
                         +" (serve zeta >= "+str(round(zeta_min_bs(wn),3))+")")
    r = np.sqrt(max(delta, 0.0))
    return (beta - r)/2, (beta + r)/2


def coefficienti_dfl(pol_1, pol_2):
    # c3, c2, c1, c0 del DFL: coefficienti di (s^2 + b1 s + a1)(s^2 + b2 s + a2), dati (b1, a1) e (b2, a2)
    p = np.polymul([1, *pol_1], [1, *pol_2])
    return p[1:]


def guadagni_bs(poli=POLI):
    # Nel formato GUADAGNI di main_hover.py: (k1, k2) per canale
    k_xy = canale_bs(*poli['xy'])
    k_att = canale_bs(*poli['att'])
    return {'x': k_xy, 'y': k_xy, 'z': canale_bs(*poli['z']),
            'phi': k_att, 'theta': k_att, 'psi': canale_bs(*poli['psi'])}


def guadagni_dfl(poli=POLI):
    # Posizione: x'''' chiuso con (canale di posizione) x (canale d'assetto); quota con la coppia in più del DFL.
    # Righe (c3, c2, c1, c0); yaw (c1, c0)
    pol = polinomi(poli)
    c_xy = coefficienti_dfl(pol['xy'], pol['att'])
    return {'x': c_xy, 'y': c_xy.copy(), 'z': coefficienti_dfl(pol['z'], pol['dfl_z']), 'psi': np.array(pol['psi'])}


def guadagni_geo(m, J, poli=POLI):
    # Guadagni diagonali fisici: posizione s^2 + (kv/m) s + kx/m, assetto s^2 + (kOmega/I) s + kR/I
    # Vanno calcolati con m e J letti dalla scena (nello scenario B cambiano)
    I = np.diag(np.array(J, dtype=float))      # solo la diagonale, come nei controllori
    pol = polinomi(poli)
    b_xy, a_xy = pol['xy']
    b_z, a_z = pol['z']
    b_att, a_att = pol['att']
    b_psi, a_psi = pol['psi']
    return {'kx': m*np.array([a_xy, a_xy, a_z]),
            'kv': m*np.array([b_xy, b_xy, b_z]),
            'kR': I*np.array([a_att, a_att, a_psi]),
            'kOmega': I*np.array([b_att, b_att, b_psi])}


def scala_poli(poli, kappa):
    # FASE 7: omega_n x kappa su tutti i canali; zeta alzato al minimo ammesso dal BS dove serve
    nuovi = {}
    for canale, (wn, zeta) in poli.items():
        wn = wn*kappa
        if canale != 'dfl_z':
            zeta = max(zeta, zeta_min_bs(wn))
        nuovi[canale] = (wn, zeta)
    return nuovi


def limiti(m, T_max):
    # Un solo dict per i tre: il BS legge angolo_max e T_max, DFL e GEO solo T_min
    return {'angolo_max': ANGOLO_MAX, 'T_max': T_max, 'T_min': FRAZIONE_T_MIN*m*G}


def parametri(nome, m, J, T_max, dt, poli=POLI, tau_f=TAU_F):
    # Argomenti del costruttore di ciascun controllore: Controllore(**parametri(nome, ...))
    comuni = {'massa': m, 'J': J, 'dt': dt, 'g': G, 'limiti': limiti(m, T_max)}
    if nome == 'BS':
        return {**comuni, 'guadagni': guadagni_bs(poli), 'tau_filtro': tau_f}
    if nome == 'DFL':
        return {**comuni, 'guadagni': guadagni_dfl(poli)}
    if nome == 'GEO':
        return {**comuni, 'guadagni': guadagni_geo(m, J, poli), 'tau_filtro': tau_f}
    raise ValueError("Controllore sconosciuto: "+str(nome)+" (BS, DFL o GEO)")


if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)

    print("\n--- Poli ---")
    for canale, (wn, zeta) in POLI.items():
        beta, alpha = polinomio(wn, zeta)
        print(canale.ljust(6)+" wn "+f"{wn:6.3f}"+"  zeta "+f"{zeta:5.3f}"+"  ->  s^2 + "+f"{beta:.4g}"+" s + "+f"{alpha:.4g}"
              +"   poli "+str(np.roots([1, beta, alpha])))

    print("\n--- BS (k1, k2) ---")
    for canale, k in guadagni_bs().items():
        print(canale.ljust(6)+str(np.array(k)))

    print("\n--- DFL (c3, c2, c1, c0); psi (c1, c0) ---")
    for canale, c in guadagni_dfl().items():
        print(canale.ljust(6)+str(c))

    print("\n--- GEO per unità di massa e d'inerzia (m = 1, J = I): moltiplicare per m e J della scena ---")
    for nome, k in guadagni_geo(1.0, np.eye(3)).items():
        print(nome.ljust(7)+str(k))

    # Test della FASE 1: i guadagni di oggi in main_hover.py
    attesi = {'x': (1.0, 1.5), 'y': (1.0, 1.5), 'z': (2.0, 3.0), 'phi': (8.0, 10.0), 'theta': (8.0, 10.0), 'psi': (4.0, 6.0)}
    ok = all(np.allclose(guadagni_bs()[c], attesi[c]) for c in attesi)
    print("\nguadagni_bs(POLI) uguale ai GUADAGNI di main_hover.py: "+("sì" if ok else "NO"))
