import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from coppeliasim_zmqremoteapi_client import RemoteAPIClient

# La classe Drone è comune a tutti gli esercizi: si trova in comune/drone.py
sys.path.append(str(Path(__file__).resolve().parents[1]))
from comune.drone import Drone
from controllore_backstepping import ControlloreBackstepping

np.set_printoptions(precision=3, suppress=True)

# --- Impostazioni della prova ---
DT = 0.005                          # passo di simulazione [s] (il default di 50 ms è troppo lento per l'assetto)
DURATA = 10.0                       # durata della prova [s]
OFFSET_HOVER = np.array([0, 0, 0])  # punto di hover rispetto alla posizione iniziale [m] (alzalo se il drone parte da terra)

# Motori (scelti da noi): con questi valori il hover cade a circa il 50% di Omega_max per una massa di 1 kg
PARAM_MOTORI = {
    'k_f': 1e-5,            # [N/(rad/s)^2]
    'k_m': 1.6e-7,          # [N*m/(rad/s)^2]
    'Omega_max': 1000.0,    # [rad/s]
    'tau_m': 0.02,          # [s]
    'segni': [1, -1, 1, -1] # eliche diagonali con lo stesso verso
}

# Guadagni (k1, k2): per guadagni grandi i poli sono circa -k1 e -k2
# assetto ~9 rad/s: almeno 5 volte sotto 1/tau_m = 50 rad/s e circa 7 volte sopra la posizione
GUADAGNI = {
    'x': (1.0, 1.5),
    'y': (1.0, 1.5),
    'z': (2.0, 3.0),
    'phi': (8.0, 10.0),
    'theta': (8.0, 10.0),
    'psi': (4.0, 6.0),
}


def controlla_drone(sim, drone):
    # Controlli una tantum prima di volare: se qualcosa qui è strano, il controllore non c'entra
    print("\n--- Controlli sul drone ---")
    print("Passo di simulazione:"+str(drone.dt)+" [s]")

    # Le terne di /base (R, omega) e di /Quadcopter (eliche, mixer) devono avere gli stessi assi
    mat = np.array(sim.getObjectMatrix(drone.drone, drone.radice)).reshape(3,4)
    print("Rotazione di base rispetto alla radice:\n"+str(mat[:,:3]))
    if not np.allclose(mat[:,:3], np.eye(3), atol=1e-3):
        print("ATTENZIONE: le due terne non sono allineate, R/omega e mixer sono in terne diverse")

    print("Posizioni delle eliche [m]:\n"+str(drone.pos_eliche))
    print("Numero di condizionamento del mixer:"+str(round(np.linalg.cond(drone.M),1)))

    Omega_hover = np.sqrt(drone.massa*9.81 / (4*drone.k_f))
    print("Omega di hover:"+str(round(Omega_hover,1))+" [rad/s] = "+str(round(100*Omega_hover/drone.Omega_max))+"% di Omega_max")
    if not 0.3 < Omega_hover/drone.Omega_max < 0.75:
        print("ATTENZIONE: rivedi k_f o Omega_max, il margine per le coppie è scarso")
    print("---------------------------\n")

    return


def grafici(log):
    t = np.array(log['t'])
    e_p = np.array(log['e_p'])
    eta = np.degrees(np.array(log['eta']))
    eta_d = np.degrees(np.array(log['eta_d']))
    tau = np.array(log['tau'])
    Omega = np.array(log['Omega'])

    fig, ax = plt.subplots(3, 2, figsize=(12, 9), sharex=True)

    ax[0,0].plot(t, e_p)
    ax[0,0].set_title('Errore di posizione [m]')
    ax[0,0].legend(['e_x', 'e_y', 'e_z'])

    for i, nome in enumerate(['phi', 'theta', 'psi']):
        ax[0,1].plot(t, eta[:,i], label=nome)
        ax[0,1].plot(t, eta_d[:,i], '--', color=ax[0,1].lines[-1].get_color(), label=nome+'_d')
    ax[0,1].set_title('Assetto [deg]')
    ax[0,1].legend()

    ax[1,0].plot(t, log['T'])
    ax[1,0].set_title('Spinta T [N]')

    ax[1,1].plot(t, tau)
    ax[1,1].set_title('Coppie in terna body [N*m]')
    ax[1,1].legend(['tau_x', 'tau_y', 'tau_z'])

    ax[2,0].plot(t, Omega)
    ax[2,0].set_title('Velocità dei motori [rad/s]')
    ax[2,0].legend(['1', '2', '3', '4'])
    ax[2,0].set_xlabel('t [s]')

    ax[2,1].plot(t, log['V_pos'], label='V posizione')
    ax[2,1].plot(t, log['V_att'], label='V assetto')
    ax[2,1].set_title('Funzioni di Lyapunov')
    ax[2,1].legend()
    ax[2,1].set_xlabel('t [s]')

    for a in ax.flat:
        a.grid(True)
    fig.tight_layout()
    plt.show()

    return


if __name__ == "__main__":
    print('Program started')

    # Connessione
    client = RemoteAPIClient()
    sim = client.require('sim')

    # Il passo si cambia a simulazione ferma, prima di creare il drone (che lo legge nel costruttore)
    sim.setFloatProperty(sim.handle_scene, 'timeStep', DT)

    drone = Drone(sim, '/Quadcopter', PARAM_MOTORI)
    controlla_drone(sim, drone)

    controllore = ControlloreBackstepping(
        massa=drone.massa,
        J=drone.J,
        guadagni=GUADAGNI,
        limiti={'angolo_max': np.radians(25), 'T_max': 0.9*drone.T_max},  # margine per le coppie
        dt=drone.dt,
    )

    # Simulazione sincronizzata: ogni sim.step() avanza di un DT
    sim.setStepping(True)
    sim.startSimulation()

    # Riferimento di hover: fermo nel punto di partenza (più l'offset), con lo yaw iniziale
    drone.aggiorna_stato()
    rif = {
        'p': drone.p + OFFSET_HOVER,
        'v': np.zeros(3),
        'a': np.zeros(3),
        'psi': drone.eta[2],
        'dpsi': 0.0,
        'ddpsi': 0.0,
    }
    print("Punto di hover:"+str(rif['p']))

    log = {k: [] for k in ['t', 'e_p', 'eta', 'eta_d', 'T', 'tau', 'Omega', 'V_pos', 'V_att']}

    try:
        while (t := sim.getSimulationTime()) < DURATA:
            # 1. stato  2. controllo  3. attuazione  4. log  5. passo
            drone.aggiorna_stato()
            T, tau = controllore.calcola(drone.p, drone.v, drone.eta, drone.omega, rif)
            drone.applica_comando(T, tau)

            log['t'].append(t)
            log['e_p'].append(controllore.e_p.copy())
            log['eta'].append(drone.eta.copy())
            log['eta_d'].append(controllore.eta_d.copy())
            log['T'].append(T)
            log['tau'].append(tau.copy())
            log['Omega'].append(drone.Omega.copy())
            log['V_pos'].append(controllore.V_pos)
            log['V_att'].append(controllore.V_att)

            if len(log['t']) % round(0.5/drone.dt) == 0:
                print("t:"+str(round(t,2))+"[s]  e_p:"+str(controllore.e_p)+"[m]  eta:"+str(np.degrees(drone.eta))+"[deg]  T:"+str(round(T,2))+"[N]")

            sim.step()
    finally:
        # La simulazione si ferma anche se qualcosa va storto (o con Ctrl+C)
        sim.stopSimulation()

    grafici(log)
