"""
Rapporto dell'analisi lineare in hover: esegue lineare.py e scrive doc/Analisi lineare in hover.html.

Il testo e i grafici stanno nel modello analisi/modello_rapporto_lineare.html; qui si calcolano i numeri
e si inseriscono nel modello al posto di /*DATI*/. Se cambiano i poli o i parametri, basta rilanciare:

    python analisi/rapporto_lineare.py
"""
import json
import subprocess
import sys
from datetime import date
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lineare as L

CARTELLA = Path(__file__).resolve().parent
MODELLO = CARTELLA / 'modello_rapporto_lineare.html'
USCITA = CARTELLA.parent / 'doc' / 'Analisi lineare in hover.html'


def r(v, cifre=4):
    # Arrotonda per il JSON: numeri finiti a 'cifre' cifre significative, infiniti e None come null
    if v is None:
        return None
    v = float(v)
    if not np.isfinite(v):
        return None
    if v == 0:
        return 0.0
    return float(f"{v:.{cifre}g}")


def lista(a, cifre=4):
    return [r(v, cifre) for v in np.asarray(a, float)]


def poli_json(p):
    return [[r(z.real), r(z.imag)] for z in np.asarray(p, complex)]


def margini_json(m):
    return {k: r(v) for k, v in m.items()}


def commit_git():
    try:
        h = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=CARTELLA, capture_output=True, text=True).stdout.strip()
        sporco = subprocess.run(['git', 'status', '--porcelain'], cwd=CARTELLA, capture_output=True, text=True).stdout.strip() != ''
        return h, sporco
    except OSError:
        return None, None


def gm_nyquist(tipo, legge, T=L.T_C, tau_m=L.TAU_M):
    # Nel discreto la fase arriva a -180° esattamente alla frequenza di Nyquist: margine di guadagno lì
    F, H, C = L.impianto_discreto(tipo, T, tau_m)
    y = C @ np.linalg.solve(-np.eye(len(F)) - F, H)
    val = -(legge.disc.K(np.array([-1.0 + 0j]))[0] @ y)
    if abs(val.imag) < 1e-9 and val.real < 0:
        return -20*np.log10(abs(val.real))
    return np.inf


def curva(nome, Lw, w=L.W, db_min=-36, db_max=42, punti=420):
    # Curva di Nichols ridotta a 'punti' campioni dentro la finestra in dB
    fase, db = L.curva_nichols(Lw, w)
    dentro = np.nonzero((db >= db_min) & (db <= db_max))[0]
    i = np.unique(np.linspace(dentro[0], dentro[-1], punti).astype(int))
    m = L.margini(Lw, w)
    return {'nome': nome, 'w': lista(w[i]), 'fase': lista(fase[i]), 'db': lista(db[i]), 'margini': margini_json(m)}


def metriche_x(t, x, th, U, rif=None):
    sotto = -np.min(x)
    return {'t_5mm': r(L.tempo_entro(t, x, 0.005)),
            'sovra_mm': r(1000*max(sotto, 0.0)),
            'incl_max': r(np.degrees(np.max(np.abs(th)))),
            'u_max': r(np.max(np.abs(U))),
            'dist_mm': None if rif is None else r(1000*np.max(np.abs(x - rif)))}


def campiona(t, *serie, ogni=4):
    return [lista(t[::ogni])] + [lista(s[::ogni]) for s in serie]


def calcola():
    D = {}
    h, sporco = commit_git()
    D['meta'] = {'data': date.today().isoformat(), 'commit': h, 'sporco': sporco,
                 'g': L.G, 'tau_m': L.TAU_M, 'T_c': L.T_C, 'tau_f': L.TAU_F,
                 'poli': {k: [r(v[0]), r(v[1])] for k, v in L.POLI.items()}}

    # --- Poli e guadagni ---
    gu = L.tabella_guadagni()
    D['guadagni'] = {c: {'wn': r(gu[c]['wn']), 'zeta': r(gu[c]['zeta']), 'beta': r(gu[c]['beta']),
                         'alpha': r(gu[c]['alpha']), 'zeta_min': r(gu[c]['zeta_min']),
                         'k': lista(gu[c]['k_bs']), 'poli': poli_json(gu[c]['poli'])}
                     for c in ['xy', 'z', 'att', 'psi']}
    D['guadagni']['dfl_xy'] = lista(gu['dfl_xy'])
    D['guadagni']['dfl_z'] = lista(gu['dfl_z'])
    D['guadagni']['dfl_psi'] = lista(gu['dfl_psi'])

    # Le tre leggi sull'asse x come vettori K: u = K (x, vx, theta, q)
    D['K_x'] = {'bs_nulle': lista(L.legge_bs_laterale(derivate='nulle').Dc),
                'bs_esatte': lista(L.legge_bs_laterale(derivate='esatte').Dc),
                'dfl': lista(L.legge_dfl_laterale().Dc)}

    # --- Leggi usate in tutto il rapporto ---
    bs_f = L.legge_bs_laterale()
    bs_n = L.legge_bs_laterale(derivate='nulle')
    bs_e = L.legge_bs_laterale(derivate='esatte')
    dfl = L.legge_dfl_laterale()
    att = L.legge_doppio(L.POLI['att'])
    quota = L.legge_doppio(L.POLI['z'])
    imb = L.legge_doppio(L.POLI['psi'])
    dfl_z = {lam: L.legge_dfl_z({**L.POLI, 'dfl_z': (lam, 1.0)}) for lam in (9.0, 15.0, 20.0, 25.0)}

    # --- Margini: continuo con mezzo passo di ritardo, e discreto esatto ---
    righe = [('Assetto φ, θ (da solo)', 'doppio', att, 'tutti'),
             ('Asse x, BS/GEO, derivata filtrata', 'laterale', bs_f, 'BS, GEO'),
             ('Asse x, BS/GEO, derivate nulle', 'laterale', bs_n, 'BS, GEO'),
             ('Asse x, DFL (= BS con derivate esatte)', 'laterale', dfl, 'DFL'),
             ('Quota z, BS/GEO', 'doppio', quota, 'BS, GEO'),
             ('Quota z, DFL, poli in più −9 doppio', 'doppio', dfl_z[9.0], 'DFL'),
             ('Quota z, DFL, poli in più −15 doppio', 'doppio', dfl_z[15.0], 'DFL'),
             ('Quota z, DFL, poli in più −25 doppio', 'doppio', dfl_z[25.0], 'DFL'),
             ('Imbardata ψ', 'doppio', imb, 'tutti')]
    D['margini'] = []
    for nome, tipo, lg, chi in righe:
        mc = L.margini(L.anello_attuatore(tipo, lg))
        md = L.margini(L.anello_attuatore_discreto(tipo, lg))
        md['gm'] = min(md['gm'], gm_nyquist(tipo, lg))
        D['margini'].append({'nome': nome, 'chi': chi, 'cont': margini_json(mc), 'disc': margini_json(md)})
    D['derivati'] = {'dpm_att': r(D['margini'][0]['disc']['pm'] - D['margini'][0]['cont']['pm']),
                     'dpm_x': r(D['margini'][1]['disc']['pm'] - D['margini'][1]['cont']['pm'])}
    D['margini_pos'] = {'nulle': margini_json(L.margini(L.anello_posizione_bs(derivate='nulle'))),
                        'filtrate': margini_json(L.margini(L.anello_posizione_bs())),
                        'filtrate_10ms': margini_json(L.margini(L.anello_posizione_bs(tau_f=0.01))),
                        'act_10ms': margini_json(L.margini(L.anello_attuatore('laterale', L.legge_bs_laterale(tau_f=0.01))))}

    # --- Carte di Nichols ---
    D['nichols'] = {
        'assetto': [curva('ωₙ = 9 (nominale)', L.anello_attuatore('doppio', att)),
                    curva('ωₙ = 10.8 (limite a 50°)', L.anello_attuatore('doppio', L.legge_doppio((L.wn_max_assetto(), 1.0)))),
                    curva('ωₙ = 14', L.anello_attuatore('doppio', L.legge_doppio((14.0, 1.0))))],
        'attuatore': [curva('BS/GEO, derivata filtrata', L.anello_attuatore('laterale', bs_f)),
                      curva('DFL', L.anello_attuatore('laterale', dfl)),
                      curva('BS/GEO, derivate nulle', L.anello_attuatore('laterale', bs_n))],
        'posizione': [curva('derivata filtrata', L.anello_posizione_bs()),
                      curva('derivate nulle', L.anello_posizione_bs(derivate='nulle'))],
        'quota': [curva('BS/GEO', L.anello_attuatore('doppio', quota)),
                  curva('DFL, −9 doppio', L.anello_attuatore('doppio', dfl_z[9.0])),
                  curva('DFL, −15 doppio', L.anello_attuatore('doppio', dfl_z[15.0])),
                  curva('DFL, −25 doppio', L.anello_attuatore('doppio', dfl_z[25.0]))],
    }

    # --- Esercizio 8: margine dell'assetto al variare di omega_n ---
    wn = np.arange(3.0, 16.01, 0.25)
    mm = [L.pm_assetto(v) for v in wn]
    D['pm_wn'] = {'wn': lista(wn), 'pm': lista([m['pm'] for m in mm]), 'gm': lista([m['gm'] for m in mm]),
                  'wc': lista([m['wc'] for m in mm]), 'wn_max': r(L.wn_max_assetto())}

    # --- Poli in anello chiuso ---
    def riga_poli(nome, tipo, lg):
        sd, raggio = L.poli_discreti(tipo, lg)
        return {'nome': nome, 'senza': poli_json(L.poli_continui(tipo, lg, tau_m=0)),
                'con': poli_json(L.poli_continui(tipo, lg)), 'disc': poli_json(sd), 'raggio': r(raggio, 6)}
    D['poli'] = [riga_poli('Asse x, BS/GEO, derivata filtrata', 'laterale', bs_f),
                 riga_poli('Asse x, BS/GEO, derivate nulle', 'laterale', bs_n),
                 riga_poli('Asse x, DFL (= BS con derivate esatte)', 'laterale', dfl),
                 riga_poli('Assetto da solo', 'doppio', att),
                 riga_poli('Quota z, BS/GEO', 'doppio', quota),
                 riga_poli('Quota z, DFL, −9 doppio', 'doppio', dfl_z[9.0]),
                 riga_poli('Imbardata ψ', 'doppio', imb)]

    # --- Gradino di 10 cm su x ---
    t, x_bd, th_bd, u_bd = L.gradino_x(bs_f, discreta=True)
    _, x_bc, th_bc, u_bc = L.gradino_x(bs_f, discreta=False)
    _, x_bn, th_bn, u_bn = L.gradino_x(bs_n, discreta=False)
    _, x_be, th_be, u_be = L.gradino_x(bs_e, discreta=False)
    _, x_dd, th_dd, u_dd = L.gradino_x(dfl, discreta=True)
    D['gradino_x'] = {
        'serie': campiona(t, 1000*x_bd, 1000*x_dd, 1000*x_bn, np.degrees(th_bd), np.degrees(th_dd), np.degrees(th_bn)),
        'metriche': [{'nome': 'BS/GEO, derivata filtrata, discreto (previsione T1-N1)', **metriche_x(t, x_bd, th_bd, u_bd, x_bd)},
                     {'nome': 'BS/GEO, derivata filtrata, continuo', **metriche_x(t, x_bc, th_bc, u_bc, x_bd)},
                     {'nome': 'BS/GEO, derivate esatte, continuo', **metriche_x(t, x_be, th_be, u_be, x_bd)},
                     {'nome': 'DFL, discreto', **metriche_x(t, x_dd, th_dd, u_dd, x_bd)},
                     {'nome': 'BS/GEO, derivate nulle, continuo', **metriche_x(t, x_bn, th_bn, u_bn, x_bd)}],
        'cont_vs_disc_mm': r(1000*np.max(np.abs(x_bd - x_bc))),
        'dfl_vs_bs_mm': r(1000*np.max(np.abs(x_bd - x_dd))),
    }

    # --- Esercizio 11: derivate esatte, filtrate, nulle ---
    der = L.confronto_derivate()
    D['derivate'] = {k: {'senza': poli_json(v['poli_senza_motore']), 'con': poli_json(v['poli_con_motore']),
                         'dist_mm': r(1000*v['dist_max']), 't_5mm': r(v['t_5mm'])} for k, v in der.items()}

    # --- Gradino di 0.5 m su z ---
    tz, z_bs, uz_bs = L.gradino_z(quota, discreta=True)
    serie_z, metr_z = [1000*z_bs], [{'nome': 'BS/GEO', 't_5': r(L.tempo_entro(tz, z_bs, 0.025)), 'dist_mm': 0.0,
                                     'sovra_mm': r(1000*max(-np.min(z_bs), 0)), 'u_max': r(np.max(np.abs(uz_bs)))}]
    for lam in (9.0, 15.0, 25.0):
        _, z, uz = L.gradino_z(dfl_z[lam], discreta=True)
        serie_z.append(1000*z)
        metr_z.append({'nome': 'DFL, −'+str(int(lam))+' doppio', 't_5': r(L.tempo_entro(tz, z, 0.025)),
                       'dist_mm': r(1000*np.max(np.abs(z - z_bs))), 'sovra_mm': r(1000*max(-np.min(z), 0)),
                       'u_max': r(np.max(np.abs(uz)))})
    D['gradino_z'] = {'serie': campiona(tz, *serie_z), 'metriche': metr_z}

    # --- Gradino di 10° sullo yaw ---
    tp, Yp, Up = L.risposta('doppio', imb, [np.radians(10), 0], 3.0)
    D['gradino_psi'] = {'t_5': r(L.tempo_entro(tp, Yp[:, 0], np.radians(0.5))), 'u_max': r(np.max(np.abs(Up))),
                        'sovra_deg': r(np.degrees(max(-np.min(Yp[:, 0]), 0)))}

    # --- Scenario C: tau_m ---
    tms = np.arange(0.005, 0.1501, 0.0025)
    serie = {'att': [], 'act': [], 'act_dfl': [], 'z_dfl': []}
    for tm in tms:
        serie['att'].append(L.margini(L.anello_attuatore('doppio', att, tau_m=tm))['pm'])
        serie['act'].append(L.margini(L.anello_attuatore('laterale', bs_f, tau_m=tm))['pm'])
        serie['act_dfl'].append(L.margini(L.anello_attuatore('laterale', dfl, tau_m=tm))['pm'])
        serie['z_dfl'].append(L.margini(L.anello_attuatore('doppio', dfl_z[9.0], tau_m=tm))['pm'])
    tab = L.margini_tau_m()
    tm_disc, tm_cont = L.tau_m_instabile()
    D['tau_m'] = {'tm': lista(1000*tms), 'serie': {k: lista(v) for k, v in serie.items()},
                  'tabella': [{'tm': r(1000*tm), **{k: margini_json(v) for k, v in m.items()}} for tm, m in tab.items()],
                  'instab_disc': r(1000*tm_disc), 'instab_cont': r(1000*tm_cont)}

    # --- Esercizio 14: kappa ---
    sk = L.sensibilita_kappa()
    D['kappa'] = []
    serie_k = []
    for k, v in sk.items():
        g = v['guadagni']
        D['kappa'].append({'kappa': k,
                           'canali': {c: {'wn': r(g[c]['wn']), 'zeta': r(g[c]['zeta']), 'beta': r(g[c]['beta']),
                                          'alpha': r(g[c]['alpha']), 'k': lista(g[c]['k_bs'])} for c in ['xy', 'z', 'att', 'psi']},
                           'dfl_xy': lista(g['dfl_xy']),
                           'pm_att': r(v['pm_att']), 'pm_act': r(v['pm_act']), 'pm_pos': r(v['pm_pos']),
                           't_5mm': r(v['t_5mm']), 'incl_max': r(v['incl_max'])})
        tk, xk, _, _ = L.gradino_x(L.legge_bs_laterale(v['poli']))
        serie_k.append(1000*xk)
    D['gradino_kappa'] = campiona(tk, *serie_k)

    # --- Scenario B e controprova a massa non nota ---
    D['massa'] = [{'k': k, 'ez_bs': r(v['ez_bs']), 'ez_dfl': r(v['ez_dfl']),
                   'att': margini_json(v['att']), 'act': margini_json(v['act'])} for k, v in L.massa_non_nota().items()]
    D['massa_dfl_25'] = [r(v['ez_dfl']) for v in L.massa_non_nota(poli={**L.POLI, 'dfl_z': (25.0, 1.0)}).values()]

    # --- Controlli contro i blocchi Verifica ---
    D['controlli'] = [{'es': es, 'descr': d, 'atteso': a, 'calcolato': r(c, 5), 'tol': tol, 'ok': bool(abs(c - a) <= tol)}
                      for es, d, a, c, tol in L.controlli()]
    return D


if __name__ == "__main__":
    dati = calcola()
    testo = MODELLO.read_text(encoding='utf-8')
    if '/*DATI*/' not in testo:
        raise SystemExit("Nel modello manca il segnaposto /*DATI*/")
    js = json.dumps(dati, ensure_ascii=False, separators=(',', ':'))
    USCITA.write_text(testo.replace('/*DATI*/', js), encoding='utf-8')

    n_ok = sum(c['ok'] for c in dati['controlli'])
    print("Rapporto scritto in "+str(USCITA))
    print("Controlli contro i blocchi Verifica: "+str(n_ok)+"/"+str(len(dati['controlli'])))
