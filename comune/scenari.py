"""
Catalogo delle prove (doc/Campagna delle simulazioni.html): 3 riferimenti x 4 scenari x 3 livelli = 36 prove
per controllore, 108 in tutto.

Una prova è un dato, non codice: descrizione(test, scenario, livello) restituisce un dict con tutto quello che
serve per eseguirla e per analizzarla (riferimento, posa iniziale, aria, carico, motori, durata, fasi del volo,
criteri di fallimento). Non contiene guadagni: quelli vengono da comune/taratura.py.

Durata, finestre e istante della raffica si calcolano dal riferimento: se cambia un livello restano giusti.

    python -m comune.scenari    -> elenca le 36 prove
"""
import numpy as np

from comune.traiettoria import Hover, Otto, TraiettoriaInterpolata

G = 9.81

# --- Impianto nominale (lo stesso di main_hover.py) ---
DT = 0.005                  # passo di simulazione e di controllo [s]
PARAM_MOTORI = {
    'k_f': 1e-5,            # [N/(rad/s)^2]
    'k_m': 1.6e-7,          # [N*m/(rad/s)^2]
    'Omega_max': 1000.0,    # [rad/s]
    'tau_m': 0.02,          # [s]
    'segni': [1, -1, 1, -1],
}

# --- Riferimenti ---
P_HOVER = [0.0, 0.0, 1.5]           # T1: punto di hover [m]
P_RILASCIO = [0.0, 0.0, 3.0]        # T1-N3: rilasciato a 3 m
CENTRO_OTTO = [0.0, 0.0, 1.5]       # T2: centro dell'otto (cx, cy, z0) [m]
A_OTTO = 1.5                        # T2: semiampiezza [m]
GIRI_OTTO = 4
WAYPOINT = [[0, 0, 1.0, 0],         # T3: quadrato di 3 m, quote alterne, yaw +90° a ogni lato (x, y, z [m], psi [rad])
            [3, 0, 2.5, np.pi/2],
            [3, 3, 1.0, np.pi],
            [0, 3, 2.5, 3*np.pi/2],
            [0, 0, 1.0, 2*np.pi]]
METODO_WAYPOINT = 'min_snap'
T_ARRIVO = 4.0                      # T3: hover nell'ultimo punto dopo la fine del riferimento [s]

# --- Livelli ---
T1_N = [{'p0': [0.1, 0.0, 0.0], 'eta0': [0, 0, 0]},                         # N1: 10 cm su x
        {'p0': [1.0, 1.0, 0.5], 'eta0': [0, 0, np.pi/2]},                   # N2: (1, 1, 0.5) m e 90° di yaw
        {'p0': [0.0, 0.0, 0.0], 'eta0': [0, np.radians(120), 0]}]           # N3: a 3 m con 120° di beccheggio
OMEGA_OTTO = [0.8, 1.4, 2.3]        # T2-N1..N3 [rad/s]: inclinazione 11.7°, 32.5°, 59.8°
V_MEDIA = [1.0, 2.0, 2.5]           # T3-N1..N3 [m/s]:   inclinazione 8.4°, 37.6°, 59.6°
LIVELLO_BASE = 2                    # A, B e C girano alla difficoltà di N2

ARIA = [(2.0, 2.0), (4.0, 4.0), (6.0, 6.0)]     # A1..A3: vento costante lungo x e raffica in più [m/s]
C_ARIA = 0.05                       # forza dell'aria per 1 m/s di velocità relativa, in unità di m g
T_RAFFICA = 2.0                     # durata della raffica 1-coseno [s]
CARICO = [1.5, 2.0, 2.5]            # B1..B3: fattore su massa e inerzia
TAU_M = [0.04, 0.07, 0.11]          # C1..C3: costante di tempo dei motori [s]

# --- Durate e fasi di T1 [s] ---
T1_DURATA = 10.0
T1_FASI = [['F1', 0.0, 6.0], ['F4', 7.0, 10.0]]
T1_A_DURATA = 16.0
T1_A_RAFFICA = 8.0
T1_A_FASI = [['F2', 4.0, 8.0], ['F3', 8.0, 13.0], ['F4', 13.0, 16.0]]

# --- Criteri di fallimento: la prova si ferma al primo che scatta ---
FALLIMENTO = {'errore_max': 3.0,    # distanza dal riferimento [m]
              'quota_min': 0.1}     # quota della terna del drone [m]: sotto, è a terra

TEST = ['T1', 'T2', 'T3']
SCENARI = ['N', 'A', 'B', 'C']
LIVELLI = [1, 2, 3]


def crea_riferimento(rif:dict):
    # Dalla descrizione del riferimento all'oggetto che lo valuta: rif(t) -> dict p, v, a, j, s, psi, dpsi, ddpsi
    if rif['tipo'] == 'hover':
        return Hover(rif['p'], rif['psi'])
    if rif['tipo'] == 'otto':
        return Otto(rif['A'], rif['omega'], rif['centro'])
    if rif['tipo'] == 'waypoint':
        return TraiettoriaInterpolata(rif['punti'], rif['metodo'], v_media=rif['v_media'])
    raise ValueError("Riferimento sconosciuto: "+str(rif['tipo']))


def descrizione(test:str, scenario:str, livello:int):
    # Tutto quello che definisce la prova (test, scenario, livello); il controllore si aggiunge al momento di eseguirla
    if test not in TEST or scenario not in SCENARI or livello not in LIVELLI:
        raise ValueError("Prova inesistente: "+str((test, scenario, livello)))

    d = {'id': test+'_'+scenario+str(livello), 'test': test, 'scenario': scenario, 'livello': livello,
         'fattore_carico': 1.0, 'tau_m': PARAM_MOTORI['tau_m'], 'aria': None,
         'fallimento': dict(FALLIMENTO)}
    l_rif = livello if scenario == 'N' else LIVELLO_BASE       # difficoltà del riferimento
    t0_raffica = None

    if test == 'T1':
        if scenario == 'N' and livello == 3:
            centro, cond = P_RILASCIO, T1_N[2]
        elif scenario == 'A':
            centro, cond = P_HOVER, {'p0': [0.0, 0.0, 0.0], 'eta0': [0, 0, 0]}   # fermo in P
        else:
            centro, cond = P_HOVER, T1_N[l_rif - 1]
        d['riferimento'] = {'tipo': 'hover', 'p': list(centro), 'psi': 0.0}
        d['posa'] = {'p': list(np.add(centro, cond['p0'])), 'eta': list(cond['eta0'])}
        if scenario == 'A':
            d['durata'], d['fasi'], t0_raffica = T1_A_DURATA, [list(f) for f in T1_A_FASI], T1_A_RAFFICA
        else:
            d['durata'], d['fasi'] = T1_DURATA, [list(f) for f in T1_FASI]

    elif test == 'T2':
        omega = OMEGA_OTTO[l_rif - 1]
        d['riferimento'] = {'tipo': 'otto', 'A': A_OTTO, 'omega': omega, 'centro': list(CENTRO_OTTO)}
        periodo = 2*np.pi/omega
        d['posa'] = {'p': list(CENTRO_OTTO), 'eta': [0, 0, 0]}          # fermo al centro, il riferimento è già in moto
        d['durata'] = GIRI_OTTO*periodo
        d['fasi'] = [['F1', 0.0, periodo]]
        for giro in range(2, GIRI_OTTO + 1):
            nome = 'F3' if (scenario == 'A' and giro == 3) else 'F2'
            d['fasi'].append([nome+'_giro'+str(giro), (giro - 1)*periodo, giro*periodo])
        t0_raffica = 4.5*np.pi/omega                                    # giro 3, al lobo x = +A

    else:   # T3
        v_media = V_MEDIA[l_rif - 1]
        d['riferimento'] = {'tipo': 'waypoint', 'punti': [list(map(float, w)) for w in WAYPOINT],
                            'metodo': METODO_WAYPOINT, 'v_media': v_media}
        tr = crea_riferimento(d['riferimento'])
        d['posa'] = {'p': list(map(float, WAYPOINT[0][:3])), 'eta': [0, 0, float(WAYPOINT[0][3])]}
        d['durata'] = tr.t_fine + T_ARRIVO
        d['fasi'] = [['F2', 0.0, tr.t_fine], ['F4', tr.t_fine, tr.t_fine + T_ARRIVO]]
        t0_raffica = float(tr.tempi[1])                                 # inizio del tratto 2, in discesa
        if scenario == 'A':
            d['fasi'].insert(1, ['F3', t0_raffica, t0_raffica + 4.0])

    if scenario == 'A':
        vento, raffica = ARIA[livello - 1]
        d['aria'] = {'c_a': C_ARIA, 'vento': [vento, 0.0, 0.0], 'raffica': [raffica, 0.0, 0.0],
                     't0': t0_raffica, 'T_g': T_RAFFICA}
    elif scenario == 'B':
        d['fattore_carico'] = CARICO[livello - 1]
    elif scenario == 'C':
        d['tau_m'] = TAU_M[livello - 1]

    return d


def prove(test=None, scenario=None, livello=None):
    # Le descrizioni della griglia, eventualmente filtrate (liste di test, scenari, livelli), in ordine di griglia
    elenco = []
    for t in TEST:
        for s in SCENARI:
            for l in LIVELLI:
                if (test is None or t in test) and (scenario is None or s in scenario) and (livello is None or l in livello):
                    elenco.append(descrizione(t, s, l))
    return elenco


def riga(d:dict):
    # Una riga leggibile per l'elenco delle prove
    rif = d['riferimento']
    if rif['tipo'] == 'hover':
        r = "hover in "+str(np.round(rif['p'], 2).tolist())+", parte da "+str(np.round(d['posa']['p'], 2).tolist()) \
            + " con eta "+str(np.round(np.degrees(d['posa']['eta']), 1).tolist())+" deg"
    elif rif['tipo'] == 'otto':
        r = "otto, omega "+str(rif['omega'])+" rad/s"
    else:
        r = "waypoint "+rif['metodo']+", v_media "+str(rif['v_media'])+" m/s"
    extra = ""
    if d['aria'] is not None:
        extra = " | vento "+str(d['aria']['vento'][0])+" m/s, raffica +"+str(d['aria']['raffica'][0])+" m/s a "+f"{d['aria']['t0']:.2f}"+" s"
    if d['fattore_carico'] != 1.0:
        extra = " | carico x"+str(d['fattore_carico'])
    if d['tau_m'] != PARAM_MOTORI['tau_m']:
        extra = " | tau_m "+str(int(round(1000*d['tau_m'])))+" ms"
    return d['id'].ljust(6)+" "+f"{d['durata']:5.2f}"+" s | "+r+extra


if __name__ == "__main__":
    elenco = prove()
    for d in elenco:
        print(riga(d))
    print(str(len(elenco))+" prove per controllore")
