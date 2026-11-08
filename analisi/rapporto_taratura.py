"""
Rapporto della taratura dei tre controllori: legge i poli da comune/taratura.py, calcola guadagni, margini e
risposte previste con analisi/lineare.py e scrive doc/Taratura dei tre controllori.html.

Tutto sul modello lineare in hover (slide 2 del PDF 09, con motore e campionamento), senza simulazioni.
Il testo sta nel modello analisi/modello_rapporto_taratura.html; qui si calcolano i numeri e si inseriscono
al posto di /*DATI*/. Se cambiano i poli in taratura.py basta rilanciare:

    python analisi/rapporto_taratura.py
"""
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np

CARTELLA = Path(__file__).resolve().parent
sys.path.insert(0, str(CARTELLA))
import lineare as L
from comune import taratura as T
from rapporto_lineare import r, lista, poli_json, margini_json, commit_git

MODELLO = CARTELLA / 'modello_rapporto_taratura.html'
USCITA = CARTELLA.parent / 'doc' / 'Taratura dei tre controllori.html'

# Guadagni di main_hover.py, scritti a mano prima della FASE 1: taratura.py deve riprodurli
GUADAGNI_MAIN_HOVER = {'x': (1.0, 1.5), 'y': (1.0, 1.5), 'z': (2.0, 3.0), 'phi': (8.0, 10.0), 'theta': (8.0, 10.0), 'psi': (4.0, 6.0)}

NOMI_CANALI = {'xy': 'x, y', 'z': 'z', 'att': 'φ, θ', 'psi': 'ψ', 'dfl_z': 'z, coppia in più del DFL'}


def margini_assi(P):
    # Margini all'attuatore, con motore e mezzo passo di ritardo (come nei blocchi Verifica)
    m = lambda tipo, legge: margini_json(L.margini(L.anello_attuatore(tipo, legge)))
    return {'att': m('doppio', L.legge_doppio(P['att'])),
            'x_bs': m('laterale', L.legge_bs_laterale(P, tau_f=T.TAU_F)),
            'x_dfl': m('laterale', L.legge_dfl_laterale(P)),
            'pos_bs': margini_json(L.margini(L.anello_posizione_bs(P, tau_f=T.TAU_F))),
            'z_bs': m('doppio', L.legge_doppio(P['z'])),
            'z_dfl': m('doppio', L.legge_dfl_z(P)),
            'psi': m('doppio', L.legge_doppio(P['psi']))}


def scelta_dfl_z(P, valori=(9.0, 15.0, 18.0, 20.0, 25.0, 30.0, 40.0)):
    # La coppia in più del DFL su z: quanto la quota si allontana da BS/GEO, margini, picco di comando
    t, z_bs, u_bs = L.gradino_z(L.legge_doppio(P['z']), durata=6.0)
    righe = []
    for lam in valori:
        Pl = {**P, 'dfl_z': (lam, 1.0)}
        legge = L.legge_dfl_z(Pl)
        _, z, u = L.gradino_z(legge, durata=6.0)
        m = L.margini(L.anello_attuatore('doppio', legge))
        righe.append({'lam': lam, 'dist_mm': r(1000*np.max(np.abs(z - z_bs))), 't_5': r(L.tempo_entro(t, z, 0.025)),
                      'pm': r(m['pm']), 'gm': r(m['gm']), 'u_max': r(np.max(np.abs(u))), 'scelta': lam == P['dfl_z'][0]})
    return {'righe': righe, 'bs': {'t_5': r(L.tempo_entro(t, z_bs, 0.025)), 'u_max': r(np.max(np.abs(u_bs)))}}


def calcola():
    P = T.POLI
    D = {}
    h, sporco = commit_git()
    D['meta'] = {'data': date.today().isoformat(), 'commit': h, 'sporco': sporco,
                 'g': T.G, 'tau_m': L.TAU_M, 'T_c': L.T_C, 'tau_f': T.TAU_F,
                 'tau_m_ms': 1000*L.TAU_M, 'T_c_ms': 1000*L.T_C, 'tau_f_ms': 1000*T.TAU_F,
                 'angolo_max': r(np.degrees(T.ANGOLO_MAX)), 'frazione_T_min': T.FRAZIONE_T_MIN}

    # --- Poli: un canale per riga ---
    D['canali'] = []
    for c in ['xy', 'z', 'att', 'psi', 'dfl_z']:
        wn, zeta = P[c]
        beta, alpha = T.polinomio(wn, zeta)
        D['canali'].append({'id': c, 'nome': NOMI_CANALI[c], 'wn': r(wn), 'zeta': r(zeta), 'beta': r(beta), 'alpha': r(alpha),
                            'zeta_min': None if c == 'dfl_z' else r(T.zeta_min_bs(wn)),
                            'poli': poli_json(L.ordina(np.roots([1, beta, alpha])))})
    D['canale'] = {c['id']: c for c in D['canali']}

    # --- Guadagni ---
    bs = T.guadagni_bs(P)
    D['bs'] = {c: lista(k) for c, k in bs.items()}
    D['bs_uguale_main'] = bool(all(np.allclose(bs[c], GUADAGNI_MAIN_HOVER[c]) for c in GUADAGNI_MAIN_HOVER))
    D['dfl'] = {c: lista(k, 6) for c, k in T.guadagni_dfl(P).items()}
    D['geo_n'] = {c: lista(k) for c, k in T.guadagni_geo(1.0, np.eye(3), P).items()}

    # Stessi poli, stessa legge: u_theta = K (x, vx, theta, q) del BS con derivate esatte e del DFL
    D['K_x'] = {'bs_esatte': lista(L.legge_bs_laterale(P, derivate='esatte').Dc),
                'dfl': lista(L.legge_dfl_laterale(P).Dc)}

    # --- Margini ---
    D['margini'] = margini_assi(P)
    D['wn_max_att'] = r(L.wn_max_assetto())
    D['scelta_dfl_z'] = scelta_dfl_z(P)

    # --- Risposte previste (modello discreto, nell'ordine del ciclo) ---
    t, x_bs, th_bs, _ = L.gradino_x(L.legge_bs_laterale(P, tau_f=T.TAU_F))
    _, x_dfl, th_dfl, _ = L.gradino_x(L.legge_dfl_laterale(P))
    _, x_es, _, _ = L.gradino_x(L.legge_bs_laterale(P, derivate='esatte'))
    ogni = 4
    D['gradino_x'] = {'t': lista(t[::ogni]), 'bs': lista(1000*x_bs[::ogni]), 'dfl': lista(1000*x_dfl[::ogni]),
                      't_5mm': r(L.tempo_entro(t, x_bs, 0.005)), 't_5mm_dfl': r(L.tempo_entro(t, x_dfl, 0.005)),
                      'sovra_mm': r(1000*max(-np.min(x_bs), 0.0)), 'incl_max': r(np.degrees(np.max(np.abs(th_bs)))),
                      'dist_mm': r(1000*np.max(np.abs(x_bs - x_dfl))), 'filtro_mm': r(1000*np.max(np.abs(x_bs - x_es)))}

    t, z_bs, _ = L.gradino_z(L.legge_doppio(P['z']), durata=4.0)
    _, z_dfl, _ = L.gradino_z(L.legge_dfl_z(P), durata=4.0)
    D['gradino_z'] = {'t': lista(t[::ogni]), 'bs': lista(1000*z_bs[::ogni]), 'dfl': lista(1000*z_dfl[::ogni]),
                      't_5': r(L.tempo_entro(t, z_bs, 0.025)), 't_5_dfl': r(L.tempo_entro(t, z_dfl, 0.025)),
                      'sovra_mm': r(1000*max(-np.min(z_bs), 0.0)), 'sovra_dfl_mm': r(1000*max(-np.min(z_dfl), 0.0)),
                      'dist_mm': r(1000*np.max(np.abs(z_bs - z_dfl)))}

    # --- Mappa dei poli: poli scelti e polo del motore ---
    mappa = []
    for c, cls in [('xy', 's1'), ('z', 's3'), ('psi', 's4'), ('att', 's2'), ('dfl_z', 'k3')]:
        p1, p2 = D['canale'][c]['poli']
        if abs(p1[1]) < 1e-6:
            # polo reale doppio: un punto solo
            mappa.append({'re': p1[0], 'im': 0.0, 'nome': NOMI_CANALI[c], 'cls': cls, 'doppio': True})
        else:
            mappa.append({'re': p1[0], 'im': abs(p1[1]), 'nome': NOMI_CANALI[c], 'cls': cls, 'doppio': False})
            mappa.append({'re': p1[0], 'im': -abs(p1[1]), 'nome': NOMI_CANALI[c], 'cls': cls, 'doppio': False, 'coniugato': True})
    mappa.append({'re': r(-1/L.TAU_M), 'im': 0.0, 'nome': 'motore, −1/τm', 'cls': 'ink', 'doppio': False, 'motore': True})
    D['mappa'] = mappa

    return D


if __name__ == "__main__":
    D = calcola()
    testo = MODELLO.read_text(encoding='utf-8')
    if '/*DATI*/' not in testo:
        raise SystemExit("Nel modello manca il segnaposto /*DATI*/")
    USCITA.write_text(testo.replace('/*DATI*/', json.dumps(D, ensure_ascii=False, separators=(',', ':'))), encoding='utf-8')
    print("Scritto "+str(USCITA))
    print("guadagni_bs uguale a main_hover.py: "+str(D['bs_uguale_main']))
    print("margini: "+json.dumps({k: (v['pm'], v['gm']) for k, v in D['margini'].items()}))
