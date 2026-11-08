"""
Esecuzione delle prove in CoppeliaSim (FASE 2 di STEPS_progetto.md).

    esegui_prova      UNA prova: (descrizione, controllore) -> log in out/<id>_<controllore>/
    esegui_campagna   le prove della griglia per un controllore, con la ripresa; la chiamano i simula.py

Il ciclo di controllo è scritto una volta sola e vale per i tre controllori: nessuna differenza tra loro può
nascere da due cicli scritti in modo diverso. I controllori non sanno in che prova sono, il drone non sa chi lo
comanda, la descrizione della prova non contiene guadagni.

Svolgimento di una prova:
  1. scena ricaricata dal file comune/drone_control.ttt: nessuna modifica della prova precedente sopravvive
  2. passo di simulazione; scenario B: massa e inerzia scalate nella scena
  3. Drone con i PARAM_MOTORI della prova (scenario C: tau_m): legge m e J, motori già a Omega_h
  4. guadagni da comune/taratura.py con m e J appena letti; controllore creato e azzerato
  5. posa iniziale della prova
  6. ciclo: stato -> riferimento -> [cronometro] calcola [cronometro] -> aria -> motori -> log -> criteri -> sim.step()
  7. simulazione ferma, log scritto tutto insieme (comune/registro.py)
"""
import argparse
import subprocess
import time
from datetime import datetime
from pathlib import Path

import numpy as np

from comune import registro, scenari, taratura
from comune.aria import Aria
from comune.drone import Drone, scala_carico
from comune.rotazioni import matrice_R

RADICE = Path(__file__).resolve().parents[1]
SCENA = Path(__file__).resolve().parent / 'drone_control.ttt'     # l'unica scena, comune ai tre controllori
PERCORSO_DRONE = '/Quadcopter'

# Grandezze interne da registrare per ciascun controllore (attributi aggiornati da calcola)
INTERNE = {
    'BS':  ['eta_d', 'mu', 'V_pos', 'V_att'],
    'DFL': ['zeta', 'v_virt', 'det_A'],
    'GEO': ['R_d', 'eta_d', 'e_R', 'e_Omega', 'Psi'],
}

# Serie temporali: le stesse chiavi in tutte le prove (doc/Campagna delle simulazioni.html, §3)
CHIAVI_SERIE = ['t', 'p', 'p_d', 'v', 'v_d', 'a_d', 'R', 'eta', 'omega', 'psi_d',
                'u_cmd', 'u_real', 'Omega_rich', 'Omega', 'saturazione', 't_calc', 'F_aria']


def versione_codice():
    try:
        h = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=RADICE, capture_output=True, text=True).stdout.strip()
        sporco = subprocess.run(['git', 'status', '--porcelain'], cwd=RADICE, capture_output=True, text=True).stdout.strip() != ''
        return {'commit': h or None, 'modifiche_non_salvate': sporco}
    except OSError:
        return {'commit': None, 'modifiche_non_salvate': None}


def versione_simulatore(sim):
    v = {}
    for nome, parametro in [('coppeliasim', 'intparam_program_full_version'), ('motore_fisico', 'intparam_dynamic_engine')]:
        try:
            v[nome] = sim.getInt32Param(getattr(sim, parametro))
        except Exception:
            v[nome] = None
    return v


def ferma(sim):
    # Ferma la simulazione e aspetta che lo sia davvero: serve prima di ricaricare la scena
    if sim.getSimulationState() != sim.simulation_stopped:
        sim.stopSimulation()
    sim.setStepping(False)
    t0 = time.time()
    while sim.getSimulationState() != sim.simulation_stopped:
        if time.time() - t0 > 10:
            raise RuntimeError("La simulazione non si ferma")
        time.sleep(0.02)

    return


def prepara_scena(sim, fattore_carico:float = 1.0):
    # Scena pulita dal file, passo di simulazione, carico dello scenario B
    ferma(sim)
    sim.loadScene(str(SCENA))
    sim.setFloatProperty(sim.handle_scene, 'timeStep', scenari.DT)
    if fattore_carico != 1.0:
        scala_carico(sim, PERCORSO_DRONE, fattore_carico)

    return


def misura_nominale(sim):
    # m e J del drone com'è nel file della scena: servono a controllare la scalatura dello scenario B
    prepara_scena(sim)
    drone = Drone(sim, PERCORSO_DRONE, scenari.PARAM_MOTORI)
    return drone.massa, drone.J.copy()


def controlla(p, p_d, T, tau, criteri:dict):
    # Primo criterio di fallimento che scatta, oppure None
    if not (np.all(np.isfinite(p)) and np.isfinite(T) and np.all(np.isfinite(tau))):
        return "valori non finiti"
    if np.linalg.norm(p_d - p) > criteri['errore_max']:
        return "errore di posizione oltre "+str(criteri['errore_max'])+" m"
    if p[2] < criteri['quota_min']:
        return "a terra: quota sotto "+str(criteri['quota_min'])+" m"
    return None


def esegui_prova(sim, descr:dict, nome:str, Classe, cartella:Path, nominale, prefisso:str = ''):
    t_reale = time.time()

    # 1-2. Scena pulita e carico
    prepara_scena(sim, descr['fattore_carico'])

    # 3. Drone con i motori della prova
    param_motori = {**scenari.PARAM_MOTORI, 'tau_m': descr['tau_m']}
    drone = Drone(sim, PERCORSO_DRONE, param_motori)

    # Se la scalatura di B o il passo non tornano, la prova non parte: sarebbe un errore silenzioso
    m0, J0 = nominale
    k = descr['fattore_carico']
    if not (np.isclose(drone.massa, k*m0, rtol=1e-6) and np.allclose(drone.J, k*J0, rtol=1e-6, atol=1e-12)):
        raise RuntimeError(descr['id']+": massa o inerzia diverse da "+str(k)+" volte le nominali ("
                           +str(drone.massa)+" kg invece di "+str(k*m0)+")")
    if abs(drone.dt - scenari.DT) > 1e-9:
        raise RuntimeError(descr['id']+": passo di simulazione "+str(drone.dt)+" s invece di "+str(scenari.DT))

    # 4. Guadagni dalla taratura, con m e J appena letti, e controllore azzerato
    parametri = taratura.parametri(nome, drone.massa, drone.J, drone.T_max, drone.dt)
    ctrl = Classe(**parametri)
    ctrl.reset()

    # 5. Riferimento, aria e posa iniziale
    rif = scenari.crea_riferimento(descr['riferimento'])
    aria = Aria(descr['aria'], drone.massa)
    drone.imposta_posa(descr['posa']['p'], matrice_R(descr['posa']['eta']))

    S = {chiave: [] for chiave in CHIAVI_SERIE}
    nomi_interne = [n for n in INTERNE[nome] if hasattr(ctrl, n)]
    I = {n: [] for n in nomi_interne}
    esito = {'riuscita': True, 'motivo': None, 't_fallimento': None}
    N = int(round(descr['durata']/drone.dt))

    # 6. Ciclo
    sim.setStepping(True)
    sim.startSimulation()
    try:
        for _ in range(N):
            t = sim.getSimulationTime()
            drone.aggiorna_stato()
            stato = {'p': drone.p.copy(), 'v': drone.v.copy(), 'R': drone.R.copy(), 'eta': drone.eta.copy(), 'omega': drone.omega.copy()}
            r = rif(t)

            # Costo di calcolo: solo la legge di controllo, senza comunicazione e simulazione
            try:
                c0 = time.perf_counter()
                T, tau = ctrl.calcola(stato, r)
                t_calc = time.perf_counter() - c0
            except np.linalg.LinAlgError as errore:
                esito = {'riuscita': False, 'motivo': "controllore singolare: "+str(errore), 't_fallimento': t}
                break

            F_aria = aria.forza(t, stato['v'])
            if aria.attiva:
                drone.applica_forza_esterna(F_aria)
            drone.applica_comando(T, tau)

            # Log: tutto all'istante t_k, prima di sim.step()
            S['t'].append(t)
            S['p'].append(stato['p'])
            S['p_d'].append(np.array(r['p'], dtype=float))
            S['v'].append(stato['v'])
            S['v_d'].append(np.array(r['v'], dtype=float))
            S['a_d'].append(np.array(r['a'], dtype=float))
            S['R'].append(stato['R'])
            S['eta'].append(stato['eta'])
            S['omega'].append(stato['omega'])
            S['psi_d'].append(r['psi'])
            S['u_cmd'].append(np.array([T, *tau], dtype=float))
            S['u_real'].append(drone.u_realizzato.copy())
            S['Omega_rich'].append(drone.Omega_richieste.copy())
            S['Omega'].append(drone.Omega.copy())
            S['saturazione'].append(drone.saturazione.copy())
            S['t_calc'].append(t_calc)
            S['F_aria'].append(F_aria)
            for n in nomi_interne:
                I[n].append(np.array(getattr(ctrl, n), dtype=float).copy())

            motivo = controlla(stato['p'], np.array(r['p']), T, tau, descr['fallimento'])
            if motivo is not None:
                esito = {'riuscita': False, 'motivo': motivo, 't_fallimento': t}
                break

            sim.step()
    finally:
        ferma(sim)

    # 7. Log
    serie = {chiave: np.array(valori) for chiave, valori in S.items()}
    interne = {n: np.array(valori) for n, valori in I.items()}
    if len(serie['t']):
        e = np.linalg.norm(serie['p_d'] - serie['p'], axis=1)
        riassunto = {'e_medio': float(np.mean(e)), 'e_max': float(np.max(e))}
    else:
        riassunto = {'e_medio': None, 'e_max': None}
    esito.update({'passi': len(serie['t']), 'passi_previsti': N, 'tempo_reale_s': time.time() - t_reale})

    metadati = {
        'id': cartella.name,
        'quaterna': {'controllore': nome, 'riferimento': descr['test'], 'scenario': descr['scenario'], 'livello': descr['livello']},
        'descrizione': descr,
        'drone': {'massa': drone.massa, 'J': drone.J, 'massa_nominale': m0, 'J_nominale': J0,
                  'param_motori': param_motori, 'T_max': drone.T_max, 'Omega_hover': drone.Omega_hover, 'DT': drone.dt},
        'controllore': {'classe': Classe.__name__, 'poli': taratura.POLI, 'tau_f': taratura.TAU_F,
                        'guadagni': parametri['guadagni'], 'limiti': parametri['limiti'],
                        'tau_filtro': parametri.get('tau_filtro')},
        'versione': {**versione_codice(), **versione_simulatore(sim), 'data': datetime.now().isoformat(timespec='seconds')},
        'esito': esito,
        'riassunto': riassunto,
        'chiavi_interne': nomi_interne,
    }
    registro.salva(cartella, metadati, serie, interne)

    testo_esito = "riuscita" if esito['riuscita'] else "FALLITA ("+esito['motivo']+" a t = "+f"{esito['t_fallimento']:.2f}"+" s)"
    testo_errore = "" if riassunto['e_max'] is None else " | e medio "+f"{1000*riassunto['e_medio']:.1f}"+" mm, e max "+f"{1000*riassunto['e_max']:.1f}"+" mm"
    print(prefisso+": "+testo_esito+testo_errore+" | "+f"{esito['tempo_reale_s']:.0f}"+" s reali")

    return metadati


def esegui_campagna(nome:str, Classe, cartella_out, argomenti = None):
    # Le prove della griglia per un controllore. argomenti: quelli della riga di comando (None = sys.argv)
    parser = argparse.ArgumentParser(description="Prove della griglia per il controllore "+nome)
    parser.add_argument('--test', nargs='+', choices=scenari.TEST, help="solo questi riferimenti")
    parser.add_argument('--scenario', nargs='+', choices=scenari.SCENARI, help="solo questi scenari")
    parser.add_argument('--livello', nargs='+', type=int, choices=scenari.LIVELLI, help="solo questi livelli")
    parser.add_argument('--rifai', action='store_true', help="rifà anche le prove già in out/")
    parser.add_argument('--elenco', action='store_true', help="elenca le prove senza simulare")
    a = parser.parse_args(argomenti)

    elenco = scenari.prove(a.test, a.scenario, a.livello)
    cartella_out = Path(cartella_out)

    if a.elenco:
        for d in elenco:
            fatta = registro.fatta(cartella_out / (d['id']+'_'+nome))
            print(scenari.riga(d)+("   [fatta]" if fatta else ""))
        print(str(len(elenco))+" prove per "+nome+", in "+str(cartella_out))
        return

    from coppeliasim_zmqremoteapi_client import RemoteAPIClient
    client = RemoteAPIClient()
    sim = client.require('sim')

    cartella_out.mkdir(exist_ok=True)
    nominale = misura_nominale(sim)
    print("Drone nominale: "+f"{nominale[0]:.4f}"+" kg. Prove da eseguire per "+nome+": "+str(len(elenco)))

    for i, d in enumerate(elenco):
        cartella = cartella_out / (d['id']+'_'+nome)
        prefisso = "["+str(i+1)+"/"+str(len(elenco))+"] "+cartella.name
        if registro.fatta(cartella) and not a.rifai:
            print(prefisso+": già fatta, la salto (--rifai per rifarla)")
            continue
        esegui_prova(sim, d, nome, Classe, cartella, nominale, prefisso)

    # Alla fine resta aperta la scena del file, senza le modifiche dell'ultima prova
    prepara_scena(sim)

    return
