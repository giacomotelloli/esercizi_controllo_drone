# Controllo Backstepping di un quadrirotore — piano di lavoro

Guida passo-passo **senza codice**: il codice lo scrivi tu.
Le formule finali sono nascoste nei blocchi **«Verifica»**: prova a ricavarle su carta prima di aprirli.

**Materiale di riferimento (già nella cartella del corso):**

- `02_Backstepping-based control/02_Backstepping.pdf` — slide Vendittelli (notazione usata qui sotto)
- `02_Backstepping-based control/BS1.pdf` — Madani & Benallegue, *Backstepping Control for a Quadrotor Helicopter* (2006)
- `02_Backstepping-based control/BS2.pdf` — Bouabdallah & Siegwart, *Full Control of a Quadrotor* (2007)
- `01_QuadModeling_2024.pdf` — modellazione

---

## 0. Python o CoppeliaSim?

**Tutti e due, in quest'ordine:**

1. **Python puro**: modello, integratore e controllore li scrivi tu. Se qualcosa non converge, sai che il problema è nel controllore o nel modello, non nel motore fisico o in una convenzione di CoppeliaSim.
2. **CoppeliaSim come "impianto"**: lo *stesso* modulo controllore Python viene collegato alla scena tramite **ZMQ Remote API**. Il drone simulato, con la sua massa, le sue inerzie e il motore fisico, diventa un test di robustezza.

Il controllore lo scrivi una volta sola: cambia solo chi fornisce lo stato e chi riceve gli ingressi.

---

## FASE 1 — Su carta (prima di scrivere qualsiasi riga)

### 1.1 Convenzioni (decidile e scrivile in cima ai tuoi appunti)

- **Frame inerziale con z verso l'alto (ENU).** È la terna $SR_I$ delle slide 5–6 ed è anche la convenzione di CoppeliaSim: usa **z-up ovunque**.
- **Angoli RPY**: $\phi$ roll, $\theta$ pitch, $\psi$ yaw, convenzione **ZYX**: $R = R_z(\psi)R_y(\theta)R_x(\phi)$.
- **Stato (12)**: $\mathbf p=(x,y,z)$, $\mathbf v=\dot{\mathbf p}$, $\eta=(\phi,\theta,\psi)$, $\omega=(p,q,r)$ velocità angolare in terna body.
- **Ingressi (4)**: spinta totale $T$ e coppie $\tau=(\tau_\phi,\tau_\theta,\tau_\psi)$.

### 1.2 Modello dinamico

**Esercizio 1** — Scrivi il modello completo in z-up: traslazione (Newton), cinematica angolare $\dot\eta = W(\eta)\,\omega$ e dinamica rotazionale (Eulero).

<details>
<summary>Verifica</summary>

$$\dot{\mathbf p} = \mathbf v,\qquad m\,\dot{\mathbf v} = -mg\,\mathbf e_3 + T\,R(\eta)\,\mathbf e_3$$

$$\dot\eta = W(\eta)\,\omega,\qquad W(\eta)=\begin{bmatrix}1 & s_\phi t_\theta & c_\phi t_\theta\\ 0 & c_\phi & -s_\phi\\ 0 & s_\phi/c_\theta & c_\phi/c_\theta\end{bmatrix}$$

$$J\,\dot\omega = \tau - \omega\times J\omega$$</details>

**Esercizio 2** — Esplicita $\ddot x,\ \ddot y,\ \ddot z$ in funzione di $\phi,\theta,\psi,T$.

<details>
<summary>Verifica</summary>

$$\ddot x = \tfrac{T}{m}(c_\phi s_\theta c_\psi + s_\phi s_\psi)$$
$$\ddot y = \tfrac{T}{m}(c_\phi s_\theta s_\psi - s_\phi c_\psi)$$
$$\ddot z = \tfrac{T}{m}\,c_\phi c_\theta - g$$

La $\ddot z$ coincide con quella della slide 5.
</details>

**Esercizio 3** — Ricava il **modello semplificato** per il progetto dell'assetto (piccoli angoli: $\dot\eta\approx\omega$, effetti giroscopici dei rotori trascurati).

<details>
<summary>Verifica</summary>

$$\ddot\phi = \dot\theta\dot\psi\,\frac{I_y-I_z}{I_x} + \frac{\tau_\phi}{I_x}\qquad
\ddot\theta = \dot\phi\dot\psi\,\frac{I_z-I_x}{I_y} + \frac{\tau_\theta}{I_y}\qquad
\ddot\psi = \dot\phi\dot\theta\,\frac{I_x-I_y}{I_z} + \frac{\tau_\psi}{I_z}$$

(cfr. eq. (7) di BS2 senza i termini aerodinamici)
</details>

> ⚠️ Il controllore si **progetta** sul modello semplificato, ma si **simula** sul modello completo (Esercizio 1). La differenza fra i due è già un primo test di robustezza.

### 1.3 La ricetta del backstepping (un canale)

Ogni canale ha la forma triangolare a doppio integratore $\dot x_1 = x_2,\ \dot x_2 = f(x) + b\,u$. I passi:

1. errore di inseguimento $e_1 = x_{1d} - x_1$
2. **controllo virtuale**: se $x_2$ fosse un ingresso, sceglierei $x_{2d} = \dot x_{1d} + k_1 e_1$, che rende $V_1=\tfrac12 e_1^2$ decrescente
3. secondo errore $e_2 = x_{2d} - x_2$
4. scrivi $\dot e_1$ e $\dot e_2$ (attenzione: in $\dot e_2$ compare $\dot x_{2d}$)
5. Lyapunov $V=\tfrac12(e_1^2+e_2^2)$: scegli $u$ in modo che $\dot V = -k_1 e_1^2 - k_2 e_2^2$
6. riscrivi $u$ in funzione di $e_1$ e $\dot e_1$ (grandezze misurabili)

**Esercizio 4** — Applica la ricetta alla **quota z** (variabile virtuale $\tilde T$, con $T=\frac{m}{c_\phi c_\theta}(g+\tilde T)$).
Verifica: devi ritrovare il risultato della **slide 6** con $K_{zp}$, $K_{zd}$.

**Esercizio 5** — Scrivi la dinamica dell'errore in anello chiuso $(e_1, e_2)$ come sistema lineare 2×2. Calcola gli autovalori. Quando sono reali? Cosa succede se $k_1=k_2$?

<details>
<summary>Verifica</summary>

$$\dot e_1 = -k_1 e_1 + e_2,\qquad \dot e_2 = -e_1 - k_2 e_2$$

$$\lambda_{1,2} = \frac{-(k_1+k_2)\pm\sqrt{(k_1-k_2)^2-4}}{2}$$

Gli autovalori sono reali se $|k_1-k_2|\ge 2$; con $k_1=k_2=k$ si ottiene $\lambda=-k\pm j$. Per guadagni grandi $\lambda\approx -k_1,\ -k_2$: le costanti di tempo sono circa $1/k_i$, utile per la taratura.
</details>

### 1.4 Posizione x, y (sottosistema sottoattuato)

$x$ e $y$ non hanno un ingresso diretto: si muovono solo inclinando la spinta. Si procede così:

1. Applica la ricetta 1.3 a $x$, $y$ e $z$ **come se avessi le accelerazioni come ingressi** → ottieni tre **ingressi virtuali** $\mu_x,\mu_y,\mu_z$ (accelerazioni desiderate).
2. **Inverti** la relazione dell'Esercizio 2: da $(\mu_x,\mu_y,\mu_z)$ e da $\psi_d$ ricava $T$, $\phi_d$, $\theta_d$.

È lo schema della slide 4: il Position Controller produce spinta e assetto desiderato, l'Attitude Controller li insegue.

**Esercizio 6** — Ricava l'inversione. *Suggerimento*: $\tfrac{T}{m}R\,\mathbf e_3 = (\mu_x,\ \mu_y,\ \mu_z+g)^T$; moltiplica a sinistra per $R_z(\psi)^T$.

<details>
<summary>Verifica</summary>

$$\mu_x = \ddot x_d + (1+k_{x1}k_{x2})e_x + (k_{x1}+k_{x2})\dot e_x\quad(\text{e analoghi per } y, z)$$

$$T = m\sqrt{\mu_x^2+\mu_y^2+(\mu_z+g)^2}$$

$$\phi_d = \arcsin\!\left(\frac{m\,(\mu_x s_\psi - \mu_y c_\psi)}{T}\right)\qquad
\theta_d = \operatorname{atan2}\!\left(\mu_x c_\psi + \mu_y s_\psi,\ \mu_z+g\right)$$
</details>

**Domande da porti (e da provare in simulazione):**

- Nell'inversione usi $\psi$ misurato o $\psi_d$? Prova entrambi.
- Per $T$: la formula "a norma" qui sopra, oppure quella della slide 6 con gli angoli **misurati** $T=\frac{m}{c_\phi c_\theta}(g+\mu_z)$? Che differenza c'è durante i transitori?
- **Saturazioni obbligatorie**: limita $|\phi_d|,|\theta_d|$ (ad es. 30°), tieni l'argomento di arcsin dentro $[-1,1]$, limita $T\in[0,\,T_{max}]$ e stai lontano da $c_\phi c_\theta \to 0$.

### 1.5 Assetto φ, θ, ψ (sottosistema completamente attuato)

**Esercizio 7** — Applica la ricetta 1.3 tre volte al modello semplificato (Esercizio 3).

<details>
<summary>Verifica</summary>

$$\tau_\phi = I_x\big[\ddot\phi_d + (1+k_{\phi1}k_{\phi2})e_\phi + (k_{\phi1}+k_{\phi2})\dot e_\phi\big] - (I_y-I_z)\dot\theta\dot\psi$$
$$\tau_\theta = I_y\big[\ddot\theta_d + (1+k_{\theta1}k_{\theta2})e_\theta + (k_{\theta1}+k_{\theta2})\dot e_\theta\big] - (I_z-I_x)\dot\phi\dot\psi$$
$$\tau_\psi = I_z\big[\ddot\psi_d + (1+k_{\psi1}k_{\psi2})e_\psi + (k_{\psi1}+k_{\psi2})\dot e_\psi\big] - (I_x-I_y)\dot\phi\dot\theta$$
</details>

**Problemi pratici da risolvere tu:**

- **Derivate degli angoli desiderati** $\dot\phi_d,\ddot\phi_d,\dot\theta_d,\ddot\theta_d$, che escono dal controllore di posizione e non dalla traiettoria. Hai tre opzioni, da provare in quest'ordine:
  1. porle a zero, il che funziona se l'anello d'assetto è molto più veloce di quello di posizione;
  2. derivarle numericamente con un filtro passa-basso ("dirty derivative");
  3. derivarle analiticamente, che è laborioso perché servono jerk e snap della traiettoria.
- **Separazione delle scale temporali**: l'anello d'assetto deve essere **almeno 4–10 volte più veloce** di quello di posizione (guadagni più alti).
- **Errore di yaw**: riportalo sempre in $[-\pi,\pi]$, altrimenti il drone fa un giro completo.
- *(Versione avanzata)* progetta sul modello completo: con $\nu$ = accelerazione angolare desiderata uscita dal backstepping, $\tau = J\,W^{-1}(\nu - \dot W\omega) + \omega\times J\omega$. Confronta con la versione semplificata.

### 1.6 Allocazione dei rotori (mixer)

**Esercizio 8** — Scrivi la matrice che lega $(\Omega_1^2,\dots,\Omega_4^2)$ a $(T,\tau_\phi,\tau_\theta,\tau_\psi)$ con il coefficiente di spinta $b$, quello di drag $d$ e il braccio $l$, poi invertila.

- Riferimenti: BS1 eq. (8) e BS2 eq. (11), entrambe in configurazione **"+"**.
- In configurazione **"x"** il braccio efficace diventa $l/\sqrt2$ e cambiano i segni.
- Saturazione: $\Omega_i\in[0,\Omega_{max}]$. Cosa succede alle coppie quando un motore satura?

Nella prima fase Python puoi saltare il mixer e applicare $T$ e $\tau$ direttamente; lo aggiungi più avanti.

### 1.7 Schema a blocchi finale (disegnalo tu)

```
traiettoria ──► [Backstepping posizione] ──μ──► [Inversione] ──T, φd, θd──► [Backstepping assetto] ──τ──► [Mixer] ──Ωi──► DRONE
 (pd, ṗd, p̈d)                                       ▲ ψd ────────────────────────────────▲                                    │
                └───────────────────── stato misurato (p, v, η, ω) ◄──────────────────────────────────────────────────────────┘
```

---

## FASE 2 — Simulazione in Python

### 2.1 Setup

- Un ambiente virtuale con `numpy` e `matplotlib` (`scipy` è opzionale).
- **Parametri**: prendi massa e inerzie dal drone della scena `drone_control.ttt` (vedi 3.2 per dove trovarli e le trappole sulle unità), con $g=9.81$. Così il modello Python e il drone di CoppeliaSim hanno gli stessi parametri e il confronto in Fase 3 è diretto.

### 2.2 Organizzazione suggerita (solo la struttura)

Separa in moduli/funzioni distinte:

| Modulo | Responsabilità |
|---|---|
| parametri | massa, inerzie, $g$, $l$, $b$, $d$, limiti, guadagni |
| modello | stato (12) + ingressi → derivata dello stato |
| traiettoria | $t$ → $\mathbf p_d,\dot{\mathbf p}_d,\ddot{\mathbf p}_d,\psi_d,\dot\psi_d,\ddot\psi_d$ |
| controllore | stato + riferimento → $T,\tau$ (posizione → inversione → assetto) |
| mixer | $T,\tau$ → $\Omega_i$ (e ritorno, con saturazioni) |
| simulazione | ciclo temporale, integrazione, logging |
| grafici | tutte le figure |

**Regola d'oro**: il controllore riceve solo lo **stato misurato** e il **riferimento** e restituisce gli **ingressi**. Non deve sapere chi c'è dall'altra parte, così lo riusi identico con CoppeliaSim.

### 2.3 Modello + test ad anello aperto (T0)

Implementa la derivata dello stato (Esercizio 1) e verifica i **segni** prima di qualsiasi controllo:

- [ ] fermo, $T=mg$, $\tau=0$ → resta fermo (gravità e spinta hanno i segni giusti)
- [ ] $T>mg$ → sale con accelerazione $(T-mg)/m$
- [ ] piccolo $\tau_\phi>0$ → $\phi$ cresce e il drone deriva verso **−y** (con $\psi=0$; verificalo con l'Esercizio 2)
- [ ] piccolo $\tau_\theta>0$ → $\theta$ cresce e il drone deriva verso **+x**
- [ ] solo $\tau_\psi$ → ruota solo lo yaw, la posizione non cambia

### 2.4 Integratore

- Usa un **passo fisso** con ingressi costanti fra due aggiornamenti del controllore (zero-order hold), come in un sistema reale. Per iniziare $dt$ = 1–5 ms.
- **RK4** scritto da te è un ottimo esercizio. In alternativa usa `solve_ivp` su ogni intervallo con ingressi costanti.
- Non chiamare il controllore *dentro* la funzione dell'ODE se ha memoria (filtri, derivate numeriche): il passo adattivo la rovina.
- Bonus: confronta con Eulero esplicito e osserva da che $dt$ in poi esplode.

### 2.5 Generatore di traiettoria

Restituisce posizione, velocità e accelerazione desiderate con **derivate analitiche**. Implementa nell'ordine:

1. **hover** (punto fisso)
2. **regolazione** come nella slide 10: partenza da $(0,0,0)$ e target $(3,3,2)$, cioè $e_x(0)=3,\ e_y(0)=3,\ e_z(0)=2$
3. **elica** come nelle slide 7–9: cerchio di raggio 1 m in $xy$, quota crescente, $\psi_d=0$
4. *(opzionale)* raccordo con polinomio di 5° grado per evitare accelerazioni iniziali enormi

### 2.6 Test incrementali (non saltarne nessuno)

| Test | Cosa attivi | Condizioni | Criterio di successo |
|---|---|---|---|
| **T1** quota | solo backstepping su z, $\tau=0$ | angoli iniziali nulli, $z_d=1$ m | $e_z\to0$; $V(t)$ **monotona decrescente** (calcolala e plottala!) |
| **T2** assetto | solo backstepping d'assetto, $T=\frac{mg}{c_\phi c_\theta}$ | $\phi_0=\theta_0=\pi/8$, $\psi_0=-\pi/8$, desiderati nulli | angoli → 0, coppie ragionevoli |
| **T3** regolazione | cascata completa | come la slide 10 | grafici simili alla colonna "backstepping" della slide 10 |
| **T4** inseguimento | cascata completa | elica | come la slide 7 / slide 11 a destra |
| **T5** disturbo | cascata completa | impulso di forza/coppia esterna a metà volo (gli scenari 2–3 delle slide sembrano mostrare proprio questo) | il drone recupera la traiettoria |
| **T6** robustezza | cascata completa | massa e inerzie nel controllore sbagliate del ±20%, rumore sulle misure, mixer con saturazione | noti un errore a regime su z con la massa sbagliata? |
| **T7** confronto | PD con compensazione di gravità | stessi scenari | riproduci il confronto delle slide 10–11 |

**Grafici da produrre sempre:**

- traiettoria 3D effettiva contro quella di riferimento
- $e_x, e_y, e_z, e_\psi$ nel tempo (come nelle slide)
- $\phi,\theta$ contro $\phi_d,\theta_d$: controlla che l'anello d'assetto insegua davvero
- $T$, $\tau$ (e $\Omega_i$ se usi il mixer), per vedere le saturazioni
- funzioni di Lyapunov $V(t)$ di ogni sottosistema

### 2.7 Taratura dei guadagni

- Usa l'Esercizio 5: per guadagni grandi i poli sono circa $-k_1$ e $-k_2$.
- Valori indicativi per partire: **posizione** $k\sim1$–$3$, **assetto** $k\sim5$–$20$, con l'assetto almeno 4 volte più veloce.
- Aumenta i guadagni finché le saturazioni non cominciano a intervenire, poi torna indietro.

---

## FASE 3 — CoppeliaSim (impianto) + controllore in Python

### 3.1 Setup

- **CoppeliaSim Edu** in una versione recente (≥ 4.6 consigliata). I nomi di alcune funzioni della Remote API cambiano fra le versioni.
- Installa il client Python: pacchetto `coppeliasim-zmqremoteapi-client` (con pip).
- Il server ZMQ (add-on) parte in automatico all'avvio di CoppeliaSim.
- Esempi ufficiali: cartella di installazione → `programming/zmqRemoteApi/clients/python`. Guardali per capire connessione e stepping.
- Apri la scena `drone_control.ttt` di questa cartella.

### 3.2 Ispezione della scena (prima di collegare il controllore)

- [ ] **Gerarchia**: individua il modello del drone e la **shape dinamica principale** a cui applicherai forze e coppie.
- [ ] **Script esistenti**: il modello *Quadcopter* standard ha uno script che insegue la sfera verde "target" con un suo PID. **Disabilitalo o eliminalo**, altrimenti combatte con il tuo controllore. Tieni pure l'oggetto target come marker visivo del riferimento.
- [ ] **Terna body**: z verso l'alto? x in avanti? L'origine coincide con il baricentro?
- [ ] **Configurazione rotori**: "+" o "x"? Dove si trova ciascuna elica e in che verso gira? Ti serve per il mixer.
- [ ] **Massa e inerzia**: somma le masse di **tutte** le shape dinamiche del modello (corpo + eliche). Leggile dal dialog *Shape dynamics properties* o via API (`sim.getShapeMass`, `sim.getShapeInertia`). ⚠️ In alcune versioni il dialog mostra le inerzie **divise per la massa**: controlla le unità. Controlla anche in quale terna è espresso il tensore d'inerzia.
- [ ] **Impostazioni di simulazione**: scegli il motore fisico (Bullet/ODE/Newton/MuJoCo) e il **passo di simulazione**. Il default di 50 ms è troppo lento per l'anello d'assetto: porta il passo a 5–10 ms. La gravità deve essere $-9.81$ su z.

### 3.3 Attuazione: tre livelli, in ordine

| Livello | Come | Pro / contro |
|---|---|---|
| **A** (consigliato per primo) | forza e coppia direttamente sul corpo: `sim.addForceAndTorque` | la cosa più vicina al tuo modello Python; niente mixer |
| **B** | una forza per ogni rotore: `sim.addForce` nella posizione di ogni elica, più una coppia di reazione per lo yaw | mette alla prova mixer e saturazioni |
| **C** | eliche originali del modello (simulate a particelle) | serve identificare la relazione comando → spinta; **non necessario** per l'esercizio |

⚠️ Per il livello A controlla nella documentazione **in che terna** vanno date forza e coppia (di norma quella assoluta). In quel caso: $F_w = T\,R\,\mathbf e_3$ e $\tau_w = R\,\tau_b$. Le forze aggiunte valgono per **un solo passo** di simulazione: vanno riapplicate a ogni ciclo.

### 3.4 Ciclo di controllo da Python (ZMQ, modalità sincrona)

I passi, da scrivere tu:

1. crea il client e ottieni l'oggetto `sim`
2. recupera gli handle degli oggetti tramite il loro path (`sim.getObject`)
3. attiva lo **stepping** (sincrono): ogni chiamata a `step` avanza la simulazione di un $dt$, quindi il comportamento è deterministico e la latenza non conta. A seconda della versione: `sim.setStepping`/`sim.step` oppure `client.setStepping`/`client.step`.
4. avvia la simulazione
5. **ciclo**: leggi il tempo di simulazione → leggi lo stato → calcola il riferimento → chiama il **tuo stesso controllore** → converti nella terna giusta → applica forza e coppia → salva i dati → `step`
6. ferma la simulazione e genera gli **stessi grafici** della Fase 2

### 3.5 Lettura dello stato: le trappole

- **Posizione**: `sim.getObjectPosition` rispetto al mondo. L'origine dell'oggetto potrebbe non essere nel baricentro.
- **Orientamento**: ⚠️ `sim.getObjectOrientation` restituisce gli angoli di Eulero di CoppeliaSim $(\alpha,\beta,\gamma)$, che **non** sono roll/pitch/yaw ZYX. Leggi la matrice (`sim.getObjectMatrix`) o il quaternione (`sim.getObjectQuaternion`) e ricava tu $R$ e $\phi,\theta,\psi$ (ottimo esercizio). In alternativa usa `sim.alphaBetaGammaToYawPitchRoll`. **Verifica** ruotando a mano il drone nella scena e stampando gli angoli.
- **Velocità**: `sim.getVelocity` sulla shape dinamica (`/Quadcopter`) restituisce velocità lineare e angolare del baricentro, calcolate dal motore fisico, **in terna mondo**. Per il tuo modello ti serve $\omega_b = R^T\omega_w$, poi $\dot\eta = W(\eta)\,\omega_b$. ⚠️ Evita `sim.getObjectVelocity`: è stimata da un passo all'altro e, senza `sim.handleflag_axis`, la parte angolare sono derivate di angoli di Eulero, non $\omega$.
- *(Opzionale)* aggiungi rumore alle misure per simulare i sensori.

### 3.6 Test in CoppeliaSim

Ripeti **T0 → T5** della Fase 2. Il primo test è applicare $F = mg$ verso l'alto e controllare che il drone stia in hover.
Sovrapponi le curve di CoppeliaSim a quelle di Python. Se differiscono, le cause tipiche sono: massa o inerzie errate, baricentro spostato rispetto all'origine, passo di simulazione, motore fisico.

### 3.7 Visualizzazione nella scena

- Sposta a ogni passo l'oggetto target lungo il riferimento (`sim.setObjectPosition`).
- Traccia la scia del drone con un drawing object (`sim.addDrawingObject`), oppure fai tutti i grafici in Python dai dati salvati.

---

## FASE 4 — Estensioni (quando tutto funziona)

- **Integral backstepping** (BS2, sezione III): aggiungi l'integrale dell'errore di posizione per annullare l'errore a regime dovuto alla massa sbagliata (cfr. T6).
- **Dinamica dei rotori** del 1° ordine nel modello Python, ad es. $G(s)=\frac{0.936}{0.178s+1}$ (BS2 eq. 6): come peggiora l'inseguimento?
- **Attrito aerodinamico** $K_t$, $K_r$ (BS1, Tabella I).
- **Confronto fra controllori** sugli stessi scenari: DFL (`03_QuadControl_DFL.pdf`) e controllo geometrico su SE(3) (cartella `05`). Vedi `09_ControllersComparison.pdf`.

---

## Checklist riassuntiva

- [ ] Fase 1 — Esercizi 1–8 risolti su carta e confrontati con le «Verifica»
- [ ] Fase 2 — T0 segni del modello verificati
- [ ] Fase 2 — T1 quota con $V(t)$ decrescente
- [ ] Fase 2 — T2 assetto
- [ ] Fase 2 — T3 regolazione (slide 10)
- [ ] Fase 2 — T4 elica (slide 7)
- [ ] Fase 2 — T5–T7 disturbi, robustezza, confronto con PD
- [ ] Fase 3 — scena ispezionata, script originale disabilitato, parametri letti
- [ ] Fase 3 — stato letto e convertito correttamente (verificato a mano)
- [ ] Fase 3 — hover con livello A, poi T1–T5 in CoppeliaSim
- [ ] Fase 3 — livello B (forze per rotore + mixer)
- [ ] Fase 4 — almeno un'estensione
