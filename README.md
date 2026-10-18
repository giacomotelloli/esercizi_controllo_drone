# esercizi_controllo_drone
esercizi di teoria del controllo per un quadrirotore 

## Backstepping con CoppeliaSim + Python — come farlo girare

Guida teorica (modello, derivazioni, test): [Backstepping_control/STEPS_backstepping.md](Backstepping_control/STEPS_backstepping.md)

La fisica (gravità, dinamica del drone) la calcola CoppeliaSim. In Python scrivi solo **traiettoria desiderata**, **controllore** e il **ciclo** che li collega alla simulazione.

### 0. Una volta sola

1. **Installa il client.** Nel terminale di VS Code: `pip install coppeliasim-zmqremoteapi-client`
2. **Apri l'esempio ufficiale:**
   `C:\Program Files\CoppeliaRobotics\CoppeliaSimEdu\programming\zmqRemoteApi\clients\python\simpleTest.py`
   Guarda la **seconda metà** (quella con `setStepping`): è esattamente lo scheletro del programma da scrivere. Il lavoro è riempire quel ciclo.
3. **Sistema la scena.** Apri `Backstepping_control/drone_control.ttt` in CoppeliaSim:
   - elimina o disabilita lo script attaccato al drone;
   - nella gerarchia a sinistra guarda come si chiama il corpo del drone (per esempio `Quadcopter`): il suo path sarà `/Quadcopter`;
   - salva.

### 1. I file (nella cartella `Backstepping_control`)

- **`traiettoria.py`**: una funzione che riceve il tempo `t` e restituisce posizione, velocità e accelerazione desiderate, più ψ desiderato. **La traiettoria si scrive qui, in Python. In CoppeliaSim non si scrive niente.** Per iniziare usa un hover: posizione (0, 0, 1), tutte le derivate a zero.
- **`controllore.py`**: una funzione (o classe) che riceve lo stato del drone e i valori desiderati, e restituisce T e τ. Qui dentro vanno le formule del backstepping.
- **`main.py`**: si collega a CoppeliaSim, avvia la simulazione, esegue il ciclo e alla fine fa i grafici.

### 2. Cosa fa `main.py`, in ordine

1. **Connessione:** crea `RemoteAPIClient()` e ottieni `sim` con `client.require('sim')`.
2. **Handle del drone:** `sim.getObject('/Quadcopter')`, con il path visto nella gerarchia.
3. **Stepping:** `sim.setStepping(True)`. Da qui CoppeliaSim avanza solo quando glielo dice lo script.
4. **Play:** `sim.startSimulation()`. **Il play lo preme lo script, non tu.**
5. **Ciclo, finché `t` è minore della durata voluta:**
   - **tempo:** `t = sim.getSimulationTime()`
   - **legge il drone** (sono gli ingressi del controllore):
     - posizione con `sim.getObjectPosition(h, sim.handle_world)`;
     - orientamento con `sim.getObjectMatrix(h, sim.handle_world)`. Restituisce 12 numeri, una matrice 3×4 letta per righe: le prime tre colonne sono R, l'ultima la posizione. Da R si ricavano φ, θ, ψ.
     - velocità con `sim.getVelocity(radice)`, dove `radice` è la shape dinamica `/Quadcopter` (non funziona su un dummy). Restituisce velocità lineare e angolare del baricentro calcolate dal motore fisico, entrambe in terna mondo; ω_body = Rᵀ·ω. ⚠️ Non usare `sim.getObjectVelocity`: è una velocità stimata da un passo all'altro e, senza `sim.handleflag_axis`, come "velocità angolare" restituisce derivate di angoli di Eulero.
   - **riferimento:** chiama la funzione di `traiettoria.py` con `t`.
   - **controllo:** chiama il controllore e ottieni T e τ.
   - **comanda il drone:**
     - converti forza e coppia: F = T · (terza colonna di R), τ_mondo = R·τ;
     - applicale con `sim.addForceAndTorque(h, F, τ_mondo)`.
   - **salvataggio:** aggiunge t, stato, desiderati, T e τ a delle liste, per i grafici.
   - **avanzamento:** `sim.step()`. CoppeliaSim calcola la fisica di un passo (gravità più la forza applicata) e il drone si muove.
6. **Stop:** `sim.stopSimulation()`. La scena torna alla posizione iniziale.
7. **Grafici** con matplotlib.

### 3. Come farlo partire

1. Apri CoppeliaSim e carica `drone_control.ttt`. **Non premere play.**
2. Nel terminale di VS Code, dentro `Backstepping_control`: `python main.py`
3. Il drone si muove nella finestra di CoppeliaSim; quando lo script finisce, si aprono i grafici.

### 4. Prima del controllore: verificare che il collegamento funzioni

- **Prova 1:** `main.py` con solo i punti 1–6, e nel ciclo solo "leggi la posizione, stampa z, `sim.step()`". Se z scende, il drone cade per gravità e **la comunicazione funziona**.
- **Prova 2:** stesso programma, ma a ogni passo applichi una forza verticale pari a massa × 9.81. Se il drone **resta fermo a mezz'aria**, lettura e comando funzionano.
  - Se la forza non ha effetto, l'handle non è quello del corpo dinamico: prendi un'altra shape del drone nella gerarchia.
- **Prova 3:** solo adesso togli la forza fissa e metti il controllore e la traiettoria.
