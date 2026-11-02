# Progetto d'esame — Confronto tra controllori per quadrirotore

Backstepping, linearizzazione in retroazione dinamica (DFL) e controllo geometrico su SE(3), confrontati in CoppeliaSim su tre riferimenti (hover, otto, waypoint) e su condizioni diverse dell'ambiente e del drone.

> **Stato: bozza.** Livelli dei disturbi, geometrie e soglie sono proposte da confermare: vedi [Decisioni](#10-decisioni).

**Materiale di riferimento** (percorsi relativi alla cartella del corso `Modulo1_Mod_Control_UAV`):

- `09_ControllersComparison.pdf` — confronto PD/BS/DFL (Vendittelli) su una spirale con attrito e raffiche, valutato a occhio sui grafici 3D. È il **contesto**, non viene replicato: qui cambiano i riferimenti, si aggiunge il controllo geometrico, si variano anche le caratteristiche del drone e tutto si misura con metriche numeriche.
- `02_Backstepping-based control-20260924/02_Backstepping.pdf` — backstepping
- `03_QuadControl_DFL.pdf`, `04_approximateDFL/` — DFL
- `05_Geometric control on SE(3)-20260924/04_QuadControl_SE(3)_2024.pdf` e `LeLoMc10.pdf` (Lee, Leok, McClamroch, *Geometric Tracking Control of a Quadrotor UAV on SE(3)*, CDC 2010) — controllo geometrico

---

## 1. Definizione del progetto

Questa sezione è la definizione formale del progetto: obiettivo, formulazione, prove, grandezze misurate e prodotti attesi. Le sezioni §2–§10 contengono i dettagli operativi.

### 1.1 Obiettivo

Confrontare in modo **quantitativo e riproducibile** tre leggi di controllo non lineari per l'inseguimento di traiettoria di un quadrirotore:

| Sigla | Metodo | Cartella | Stato |
|---|---|---|---|
| **BS** | Backstepping su angoli di Eulero | `Backstepping_control/` | scritto, provato in hover |
| **DFL** | Linearizzazione in retroazione dinamica (compensatore su $T, \dot T$) | `DFL_control/` | da scrivere |
| **GEO** | Controllo geometrico su SE(3) | `Geometric_control/` | da scrivere |

per stabilire **in quali condizioni le loro prestazioni differiscono** e ricondurre ogni differenza osservata alla struttura del metodo.

Il confronto avviene a parità di impianto, riferimento, condizioni iniziali e taratura, al variare di due fattori: la **difficoltà del riferimento** e le **condizioni operative**, cioè l'ambiente e le caratteristiche fisiche del drone.

### 1.2 Formulazione

**Impianto**, uguale per tutti i controllori. La dinamica del corpo rigido la calcola CoppeliaSim; motori e mixer sono in `comune/drone.py`:

$$
\begin{aligned}
\dot p &= v, & m\,\dot v &= -mg\,e_3 + T\,R\,e_3 + F_{aria} \\
\dot R &= R\,\hat\omega, & J\,\dot\omega &= -\omega \times J\omega + \tau \\
\begin{pmatrix} T \\ \tau \end{pmatrix} &= M\,\Omega^2, & \tau_m\,\dot\Omega_i &= \Omega_{cmd,i} - \Omega_i, \qquad 0 \le \Omega_i \le \Omega_{max}
\end{aligned}
$$

con $\Omega^2$ inteso componente per componente e $M$ la matrice di allocazione delle quattro eliche.

**Uscite e riferimento.** Si inseguono le uscite piatte $y = (p, \psi)$. Il riferimento $y_d(t) = (p_d(t), \psi_d(t))$ è fornito con le derivate fino a $p_d^{(4)}$ e $\ddot\psi_d$.

**Controllori.** Ciascuno è una legge

$$
u = (T_{cmd}, \tau_{cmd}) = \kappa_c\big(x,\ \xi_c,\ y_d(t), \dots, y_d^{(4)}(t)\big), \qquad c \in \{\text{BS}, \text{DFL}, \text{GEO}\}
$$

dove $x = (p, v, R, \omega)$ è lo stato misurato e $\xi_c$ l'eventuale stato interno del controllore, cioè ciò che il controllore si porta dietro da un passo all'altro:

- **BS**: il filtro che stima le derivate di $\phi_d$ e $\theta_d$;
- **DFL**: il compensatore dinamico $\zeta = (T, \dot T)$, che fa parte del metodo;
- **GEO**: nessuno stato interno, se $\Omega_d$ e $\dot\Omega_d$ si calcolano analiticamente da jerk e snap. In quel caso è una legge statica: il comando dipende solo dallo stato e dal riferimento nell'istante corrente. Se invece si ricavano derivando numericamente $R_d$, serve un filtro come nel BS.

Nessuno dei tre modella la dinamica dei motori, le saturazioni o $F_{aria}$.

**Prova.** Una prova è individuata dalla quaterna $(c, r, s, \ell)$:

- $c$ controllore;
- $r \in \{T1, T2, T3\}$ riferimento;
- $s \in \{N, A, B, C\}$ scenario;
- $\ell$ livello: l'intensità dello scenario, oppure, nel nominale, la difficoltà del riferimento (condizione iniziale, $\omega$, `v_media`).

Per ogni fase del volo $f \in \{F1, F2, F3, F4\}$, la prova produce un vettore di metriche $m(c, r, s, \ell, f)$. Su questi vettori si fanno due confronti:

- **Confronto tra metodi.** A $(r, s, \ell, f)$ fissati cambia solo $c$: ogni differenza tra i tre vettori è attribuibile al controllore.
- **Effetto di uno scenario.** A $c$ fissato, confrontare $m(c, r, s, \ell, f)$ con il nominale isola l'effetto della sola condizione cambiata.

### 1.3 Domande e ipotesi

1. A bassa difficoltà, a parità di poli, i tre metodi sono equivalenti?
2. Come cambia la classifica quando il riferimento diventa più difficile: più veloce, con angoli più grandi, meno regolare?
3. Chi recupera da condizioni iniziali estreme (drone inclinato o capovolto)?
4. Chi resiste meglio a vento e raffiche?
5. Chi se la cava meglio con un drone che ha meno riserva di spinta, e chi soffre di più con motori lenti?
6. Quanto costa calcolare ciascuna legge?

Le ipotesi sono scritte prima di simulare, per non leggere i risultati a posteriori: alla fine ognuna sarà **confermata** o **smentita** con i numeri.

| # | Ipotesi | Perché ce la si aspetta | Dove si verifica |
|---|---|---|---|
| H1 | A bassa difficoltà, a parità di poli, le differenze sono piccole (RMSE entro il 10%) | vicino all'hover le tre leggi linearizzate coincidono; nelle slide 09 BS e DFL, senza disturbi, sono quasi sovrapposti | T1-N1, T2 a ω = 0.8, T3 a 1 m/s |
| H2 | Il DFL è il più sensibile ai disturbi non modellati | cancella esattamente il modello nominale; nelle slide 09 con le raffiche la sua deviazione appare maggiore | colonna A |
| H3 | GEO recupera da qualsiasi inclinazione iniziale; BS e DFL falliscono quando il recupero passa per $\theta = \pm 90°$ (beccheggio), non sul rollio | GEO non usa coordinate locali; con ZYX la singolarità di $W(\eta)$ è solo sul beccheggio | T1-N3 |
| H4 | Vicino al limite di spinta il DFL degrada per primo | nel DFL la spinta è uno stato del compensatore ($\zeta$), non una misura: quando i motori saturano, $\zeta$ si stacca dalla spinta vera, come un integratore in windup. BS e GEO sono leggi statiche | colonna B |
| H5 | Con motori lenti degradano tutti a livelli simili; se il DFL cede prima è per la sua ipotesi di spinta istantanea | i poli d'assetto sono gli stessi per tutti (§2.1); il DFL usa $\ddot T$ come ingresso | colonna C |
| H6 | Con un riferimento poco regolare il BS degrada meno degli altri | il BS usa solo l'accelerazione; GEO ha bisogno del jerk, il DFL dello snap | T3, variante di regolarità |
| H7 | Costo di calcolo BS < GEO < DFL, tutti molto sotto il passo | il DFL inverte una matrice 4×4 piena di termini trigonometrici | tutte le prove |

### 1.4 Le prove

**Fattore 1: il riferimento.** Tre riferimenti, ciascuno con un parametro che ne regola la difficoltà (dettagli in §5.4–§5.6):

| | Riferimento | Parametro di difficoltà | Intervallo |
|---|---|---|---|
| **T1** hover | $p_d$ costante | condizione iniziale: errore di posizione, inclinazione iniziale su rollio e beccheggio | fino a 1.5 m e 179° |
| **T2** otto | lemniscata di Gerono, $A = 1.5$ m | pulsazione $\omega$ | da 12° a 47° di inclinazione richiesta |
| **T3** waypoint | percorso 3D `min_snap` con spigoli, salite, discese e rotazioni di yaw | `v_media`, più una variante di regolarità (`cubica`, `min_jerk`, `min_snap`) | da 8° a 60° |

**Fattore 2: lo scenario.** Una condizione alla volta, con tre livelli uguali per tutti i riferimenti (dettagli in §5.3):

| | Scenario | Cosa cambia | Livelli |
|---|---|---|---|
| **N** | nominale | niente | — |
| **A** | aria | vento costante + raffica, con attrito sulla velocità relativa all'aria | forza costante da 0.1 a 0.3 $mg$, di picco da 0.2 a 0.6 $mg$ |
| **B** | carico | massa e inerzia, note al controllore; `Omega_max` resta fisso | ×1.5, ×2, ×2.5, cioè $T/W$ = 2.67, 2.0, 1.6 |
| **C** | motori lenti | costante di tempo dei motori | $\tau_m$ = 40, 60, 80 ms (nominale 20 ms) |

Nel nominale si scandisce la difficoltà del riferimento. Gli scenari A, B, C girano a una difficoltà fissa: per T1 il gradino grande, per T2 e T3 il livello più alto a cui, in nominale, tutti e tre i controllori riescono con margine. In tutto sono circa 49 prove per controllore (§5.8).

**Parità** (dettagli in §2): stesso impianto, passo di 5 ms, stesso stato misurato, stesso riferimento, stesse condizioni iniziali (drone fermo, motori già alla velocità di hover), solo le saturazioni fisiche dei motori. **Taratura** per assegnazione dei poli: stessi poli nominali del sistema linearizzato attorno all'hover, verificati con un gradino di 10 cm (T1-N1) su cui le tre risposte devono quasi sovrapporsi.

### 1.5 Grandezze misurate

Si misurano solo le uscite comuni ai tre controllori, posizione e yaw, separatamente per fase del volo: F1 convergenza, F2 regime o inseguimento, F3 disturbo, F4 arrivo (§3). Le definizioni esatte sono in §4.

| Categoria | Grandezze | A quale domanda risponde |
|---|---|---|
| Accuratezza | errore medio, RMSE, errore massimo (totale e per asse), errore di yaw, errore a regime | quanto bene insegue |
| Tempi | tempo di convergenza $t_c$, sovraelongazione, picco e tempo di recupero $t_r$ dopo un disturbo, ritardo d'arrivo, quota persa | quanto in fretta converge e recupera |
| Sforzo | RMS di $T - mg$ e di $\|\tau\|$, regolarità dei comandi, % di passi in saturazione, inclinazione massima | quanto spende e quanto vicino ai limiti lavora |
| Costo di calcolo | tempo di una chiamata al controllore: mediana, 99° percentile, massimo, rapporto con il passo | se gira in tempo reale, e con quale margine |
| Esito | successo o fallimento; nelle scansioni, la soglia di rottura | fin dove funziona |

**Criterio di rilevanza**, fissato prima di vedere i risultati: una differenza conta se supera il 10% del valore migliore, 5 mm in assoluto e la variabilità tra prove ripetute o tra giri (§6.4).

### 1.6 Prodotti attesi

- Tabelle per riferimento, scenario, livello e fase (§6.3).
- Una mappa di sintesi 3×4 (riferimenti × scenari) con il controllore migliore in ogni cella e il distacco dal secondo.
- Curve di difficoltà (errore in funzione dell'inclinazione richiesta) e curve di degrado (errore in funzione del livello dello scenario), con la soglia di rottura di ogni controllore.
- Per ciascuna ipotesi H1–H7, un verdetto motivato con un numero o un grafico.
- Codice e log riproducibili: ogni log contiene i parametri della prova e l'hash del commit.

---

## 2. Condizioni di parità

Una differenza tra i controllori significa qualcosa solo se tutto il resto è uguale.

| Cosa | Scelta | Note |
|---|---|---|
| Impianto | stessa scena `.ttt`, stessa classe `Drone`, stessi `PARAM_MOTORI` (salvo lo scenario C) | la dinamica dei motori non è modellata da nessuno dei tre: è un errore di modello comune |
| Passo | `DT = 0.005` s per fisica e controllo | |
| Stato misurato | uguale per tutti: $p, v, R, \eta, \omega$ da `Drone.aggiorna_stato()` | GEO usa $R$, BS e DFL usano $\eta$, ricavati dalla stessa matrice |
| Riferimento | stesso oggetto traiettoria: $p, v, a, j, s, \psi, \dot\psi, \ddot\psi$ | il DFL serve fino allo snap, GEO usa jerk e snap per $\Omega_d, \dot\Omega_d$ |
| Condizioni iniziali | il drone viene messo dallo script nella posa iniziale del test, fermo, con i **motori già alla velocità di hover** | oggi `Drone.Omega` parte da zero: ogni prova comincerebbe con una caduta di qualche centimetro, più lunga con i motori lenti di C. Per il DFL: $\zeta(0) = (mg, 0)$ |
| Saturazioni | solo quelle fisiche dei motori (`0 ≤ Ω ≤ Omega_max`), identiche | vedi nota sui limiti del BS |
| Parametri di modello | il controllore riceve $m$ e $J$ letti dalla scena, quindi in B **conosce** la massa nuova | B misura l'agilità del drone, non l'errore di modello; l'errore di modello è una controprova (§9) |
| Taratura | stessi poli nominali | §2.1 |

**Limiti del BS.** Oggi il BS limita l'inclinazione a 25° (`angolo_max`), la spinta al 90% di `T_max` e impone $F_z \ge 0.1\,mg$. Gli altri due non hanno questi limiti. O si applicano gli stessi limiti a tutti, o li si rende non attivi (per esempio `angolo_max` = 60°). Altrimenti su T2 e T3, dove servono fino a 60°, il BS perde per un limite che gli altri non hanno, non per il metodo.

**Derivate dell'assetto desiderato nel BS.** Con `tau_filtro=None` le derivate di $\phi_d, \theta_d$ sono poste a zero, quindi il BS insegue l'assetto in ritardo. DFL e GEO invece usano jerk e snap esatti. Su T2 e T3, dove l'assetto cambia in fretta, si misurerebbe questa approssimazione e non il metodo (nelle slide il BS ha un termine di *feedforward*). **Scelta proposta:** come base la derivata filtrata, che è già nella classe (`tau_filtro`), da dichiarare nella relazione. La versione con derivate analitiche, ricavate da jerk e snap, è un miglioramento (§9).

**Convenzioni.** Il repo usa z verso l'alto. Le slide 09 le mescolano (z verso il basso per PD e DFL, verso l'alto per il BS) e `LeLoMc10.pdf` ha $b_3$ verso il basso: quando si trascrivono le formule, controllare i segni di $g$ e di $T$.

### 2.1 Taratura: stessi poli

Vicino all'hover ogni asse laterale è una catena del quart'ordine: posizione → velocità → inclinazione → velocità angolare. Ciascun controllore piazza questi quattro poli a modo suo, ma si possono **mettere nello stesso posto**:

- **BS**: il canale $(k_1, k_2)$ dà l'errore $\ddot e + (k_1 + k_2)\dot e + (1 + k_1 k_2)e = 0$. La posizione dà la coppia di poli lenta, l'assetto quella veloce.
- **DFL**: per $x, y, z$ il polinomio d'errore è del quart'ordine, con coefficienti $c_0 \dots c_3$. Si scelgono come i coefficienti di $(s - \lambda_1)(s - \lambda_2)(s - \lambda_3)(s - \lambda_4)$, con $\lambda_{1,2}$ = poli di posizione del BS e $\lambda_{3,4}$ = poli d'assetto del BS. Anche $z$ è del quart'ordine (il compensatore rende $\ddot T$ l'ingresso): i due poli in più vanno alla velocità dell'assetto. Per $\psi$ (second'ordine) si usano i poli di yaw del BS.
- **GEO**: vicino all'hover la posizione dà $s^2 + \frac{k_v}{m}s + \frac{k_x}{m}$ e l'assetto, asse per asse, $s^2 + \frac{k_\Omega}{J_i}s + \frac{k_R}{J_i}$. Con guadagni diagonali ($K_x$, $K_v$, $K_R$, $K_\Omega$) si ottengono esattamente gli stessi poli del BS, anche su $z$ e sullo yaw.

**Punto di partenza: i guadagni attuali del BS** (`main_hover.py`):

| Canale | $(k_1, k_2)$ | Polinomio | Poli |
|---|---|---|---|
| $x, y$ | (1.0, 1.5) | $s^2 + 2.5s + 2.5$ | $-1.25 \pm 0.97j$ ($\omega_n$ ≈ 1.6 rad/s) |
| $z$ | (2.0, 3.0) | $s^2 + 5s + 7$ | $-2.5 \pm 0.87j$ |
| $\phi, \theta$ | (8, 10) | $s^2 + 18s + 81$ | $-9$ doppio |
| $\psi$ | (4, 6) | $s^2 + 10s + 25$ | $-5$ doppio |

L'equivalenza è esatta per il DFL e approssimata per BS e GEO (vale finché assetto e posizione hanno tempi ben separati; qui il rapporto è circa 6). Per questo serve un **test di taratura**, T1-N1: un gradino di 10 cm in hover. Le tre risposte devono quasi sovrapporsi; se non succede, la taratura non è equa e va corretta prima di qualsiasi altro confronto.

---

## 3. Fasi del volo

Ogni test dichiara **in anticipo** le proprie finestre temporali, uguali per tutti i controllori: le metriche si calcolano sugli stessi intervalli. Il tempo di convergenza è una misura *dentro* la fase, non ne definisce i confini.

| Fase | T1 hover | T2 otto | T3 waypoint |
|---|---|---|---|
| **F1** convergenza | dall'inizio alla fine del transitorio | primo giro (il drone parte fermo, il riferimento no) | — (parte fermo sul primo punto, senza errore) |
| **F2** regime / inseguimento | in A: sotto vento costante, prima della raffica | giri 2–4, metriche **giro per giro** | ogni tratto tra due waypoint |
| **F3** disturbo | in A: raffica + recupero | in A: il giro con la raffica, confrontato con quello prima | in A: il tratto con la raffica |
| **F4** arrivo / regime finale | ultimi secondi | — | dalla fine della traiettoria alla fine della prova |
| **TOT** | intera prova | intera prova | intera prova |

Le finestre in secondi sono in §5, test per test.

---

## 4. Metriche

Notazione: $e_p(t) = p_d(t) - p(t)$, $\|e_p\|$ norma euclidea, $e_\psi = \mathrm{wrap}(\psi_d - \psi)$ in $[-\pi, \pi]$, campioni $t_k$ con $k = 1 \dots N$ dentro la fase.

**Si confrontano solo le uscite comuni**: posizione e yaw, cioè le uscite piatte $(x, y, z, \psi)$ che tutti e tre inseguono. L'errore d'assetto rispetto all'assetto desiderato **non** è confrontabile, perché ogni controllore si calcola il proprio ($\eta_d$ per BS e DFL, $R_d$ per GEO): va nel log solo per la diagnosi.

### 4.1 Accuratezza

| Metrica | Definizione | Unità |
|---|---|---|
| Errore medio | $\bar e = \frac{1}{N}\sum_k \|e_p(t_k)\|$ | m |
| RMSE | $\sqrt{\frac{1}{N}\sum_k \|e_p(t_k)\|^2}$ | m |
| Errore massimo | $\max_k \|e_p(t_k)\|$ | m |
| Errore massimo per asse | $\max_k \lvert e_{p,i}(t_k)\rvert$, $i = x, y, z$ (mostra su quale asse agisce un disturbo) | m |
| Errore di yaw | RMS e massimo di $\lvert e_\psi \rvert$ | deg |
| Errore a regime | media di $\|e_p\|$ nella finestra di regime | m |
| IAE | $\int \|e_p\|\,dt$ (combina ampiezza e durata) | m·s |
| Errore ai waypoint (T3) | $\|e_p\|$ negli istanti di passaggio previsti dai waypoint | m |
| Errore di contorno (facoltativo, T2 e T3) | distanza di $p(t)$ dal punto più vicino del percorso: distingue "segue la forma ma in ritardo" da "taglia le curve" | m |

Su T2 le metriche di F2 si calcolano **per ogni giro** e si riportano come media ± deviazione standard sui giri 2–4: si ottiene una dispersione anche da simulazioni deterministiche.

### 4.2 Tempi di convergenza

- **Tempo di convergenza** $t_c$: il primo istante *dopo il quale* $\|e_p\|$ resta **sempre** sotto la soglia $\varepsilon$ fino alla fine della fase, contato dall'inizio della fase. Non il primo attraversamento: un controllore che oscilla attorno alla soglia sembrerebbe veloce. Se l'errore non resta mai sotto soglia: "non converge".
- **Soglia** $\varepsilon$: 5 cm in posizione, 2° in yaw. $t_c$ dipende molto da $\varepsilon$: riportarlo anche per 2 cm e 10 cm.
- **Sovraelongazione** (solo T1): con $d = e_p(0)$, il gradino iniziale, e $\hat d = d / \|d\|$,
  $M_p = \max\left(0,\ \max_k \frac{-e_p(t_k) \cdot \hat d}{\|d\|}\right)$, in percentuale. È la parte del gradino di cui il drone "passa oltre" il bersaglio.
- **Quota persa** (T1-N3): $z(0) - \min_k z(t_k)$, quanto scende il drone mentre si raddrizza.
- **Picco di deviazione** (F3): $\max_k \|e_p\|$ dentro la finestra del disturbo.
- **Tempo di recupero** $t_r$ (F3): dall'inizio del disturbo al primo istante dopo il quale $\|e_p\| < \varepsilon$ per sempre.
- **Ritardo d'arrivo** (F4): $t_{arrivo} - t_{fine}$, dove $t_{arrivo}$ è il primo istante dopo il quale il drone resta entro $\varepsilon$ dal punto finale.

### 4.3 Sforzo di controllo e attuazione

| Metrica | Definizione |
|---|---|
| Spinta | RMS di $T - mg$ |
| Coppie | RMS di $\|\tau\|$ |
| Regolarità | RMS di $\Delta T / \Delta t$ e di $\|\Delta\tau\| / \Delta t$: mette in evidenza picchi e chattering |
| Saturazione | % di passi con almeno un motore saturo (a 0 o a `Omega_max`) |
| Inclinazione massima | $\max_k \arccos(R_{33})$, angolo tra l'asse del drone e la verticale [deg] |

Servono due modifiche a `Drone`: esporre un **flag di saturazione** (oggi il taglio avviene dentro `applica_comando` e `comanda_motori` senza lasciare traccia) e permettere di loggare $T, \tau$ **comandati** (dal controllore) e **realizzati** (da $\Omega$ tramite la matrice $M$), che in saturazione sono diversi.

### 4.4 Costo di calcolo ("tempo di esecuzione")

- Si misura **solo la chiamata al controllore** (`calcola`) con `time.perf_counter()`. La comunicazione ZMQ e `sim.step()` sono esclusi: dominano il tempo totale, sono uguali per tutti e nasconderebbero le differenze. Per lo stesso motivo il tempo totale della simulazione **non** è una metrica.
- Si riportano mediana, 99° percentile e massimo per passo, in µs, e il rapporto con il passo (`t_calc / DT`).
- In Python le misure sono rumorose: niente `print` nel ciclo, stessa macchina, più ripetizioni. Più affidabile: un **benchmark offline**, che ripassa gli stati salvati nel log attraverso il controllore, fuori da CoppeliaSim e molte volte. Così tutti e tre ricevono gli stessi ingressi, senza interferenze del simulatore.
- Limite da dichiarare: si misura l'implementazione Python/numpy, non l'algoritmo. Con vettori così piccoli pesa soprattutto il costo fisso di ogni chiamata numpy.

### 4.5 Esito della prova

- Una prova **fallisce** se $\|e_p\| > 2$ m, se l'inclinazione supera 80° fuori da T1-N3, se il drone tocca terra, o se compaiono NaN. La prova si ferma; le metriche delle fasi non raggiunte valgono "—".
- Per le scansioni (velocità, intensità del disturbo, inclinazione iniziale) l'esito diventa una **soglia di rottura**: il primo livello a cui il controllore fallisce.

### 4.6 Quali metriche in quale fase

| Metrica | F1 | F2 | F3 | F4 | TOT |
|---|:-:|:-:|:-:|:-:|:-:|
| Errore medio, RMSE, massimo | ✓ | ✓ | ✓ | ✓ | ✓ |
| Massimo per asse | | ✓ | ✓ | | |
| Errore di yaw | ✓ | ✓ | ✓ | ✓ | |
| Errore ai waypoint, di contorno | | ✓ | | | |
| $t_c$, sovraelongazione, quota persa | ✓ | | | | |
| Picco, $t_r$, IAE | | | ✓ | | |
| Errore a regime, ritardo d'arrivo | | ✓ | | ✓ | |
| Sforzo, saturazione, inclinazione massima | ✓ | ✓ | ✓ | | ✓ |
| Costo di calcolo, esito | | | | | ✓ |

---

## 5. Test e scenari

### 5.1 Lo schema

|  | **N** nominale | **A** aria | **B** carico | **C** motori lenti |
|---|---|---|---|---|
| **T1 hover** | condizioni iniziali: gradino piccolo, gradino grande, drone inclinato | resistere a vento e raffica | gradino grande con drone pesante | gradino grande con motori lenti |
| **T2 otto** | scansione della velocità | a velocità di riferimento | a velocità di riferimento | a velocità di riferimento |
| **T3 waypoint** | scansione della velocità + variante di regolarità | a velocità di riferimento | a velocità di riferimento | a velocità di riferimento |

- **N** è il termine di paragone di ogni riga: senza, l'effetto di uno scenario non si può attribuire.
- **A, B, C** cambiano una sola condizione, con **gli stessi tre livelli in tutti i test** (lieve, medio, forte). Così le colonne sono confrontabili: si vede come lo stesso disturbo pesa su un hover e su un otto veloce.
- Su T2 e T3 gli scenari girano a una **difficoltà di riferimento**, scelta *dopo* la scansione nominale: il livello più alto a cui tutti e tre riescono con margine. Così l'effetto dello scenario non si confonde con un fallimento che c'era già in nominale.

### 5.2 Perché massa e motori

L'agilità del drone, cioè la capacità di inclinarsi, ruotare e reagire di colpo, dipende da due limiti diversi:

- **Quanto può spingere**: il rapporto spinta/peso $T_{max}/(mg)$. Con i parametri attuali l'hover è al 50% di `Omega_max`, quindi $T_{max} = 4\,mg$. Fissa l'accelerazione massima, e quindi l'inclinazione sostenibile e la riserva per correggere.
- **Quanto in fretta può cambiare la spinta**: la costante di tempo dei motori $\tau_m$. Fissa la rapidità con cui si cambia assetto, cioè i comportamenti impulsivi.

Massa e `Omega_max` agiscono **sulla stessa grandezza**. Se $J$ cresce in proporzione a $m$, un drone $k$ volte più pesante si comporta, in accelerazioni, *esattamente* come quello originale con $T_{max}/k$. Basta quindi variarne uno: B usa la massa, C usa $\tau_m$. Così le due colonne "drone" coprono i due limiti indipendenti.

Una conseguenza meno ovvia. A spinta totale costante, la coppia si ottiene alzando i motori da un lato di $\delta$ (in $\Omega^2$) e abbassando l'altro della stessa quantità, quindi $\delta_{max} = \min(\Omega_h^2,\ \Omega_{max}^2 - \Omega_h^2)$. Finché l'hover sta sotto il 71% di `Omega_max` ($T/W > 2$), il limite è il motore che va **a zero**. Con $J \propto m$, fino a massa ×2 l'accelerazione angolare massima **non cambia**: il carico toglie solo riserva di spinta. Oltre, il limite diventa il motore che va **al massimo**, e a ×2.5 l'accelerazione angolare massima scende del 40%. Per questo il livello forte di B sta oltre ×2.

### 5.3 Le perturbazioni (uguali per tutti i test)

**A — Aria.** Un unico modello per attrito, vento e raffiche: una forza applicata a ogni passo sul corpo (`radice`), nel baricentro, con `sim.addForceAndTorque`:

$F_{aria} = -c_a\,\big(v - v_{vento}(t)\big)$

con $c_a$ tale che 1 m/s di velocità relativa dia 0.05 $mg$. Il vento è costante lungo X; la raffica è un profilo "1-coseno" che si somma al vento: $\Delta v(t) = \frac{\Delta v_{max}}{2}\left(1 - \cos\frac{2\pi(t - t_0)}{T_g}\right)$ per $t_0 \le t \le t_0 + T_g$, con $T_g = 2$ s. In nominale non c'è aria ($c_a = 0$). Su T2 e T3 il drone sente attrito anche per il proprio moto (a 3 m/s: 0.15 $mg$).

| Livello | Vento | Raffica (in più) | Forza a drone fermo: costante / di picco |
|---|---|---|---|
| A1 lieve | 2 m/s | +2 m/s | 0.1 / 0.2 $mg$ |
| A2 medio | 4 m/s | +4 m/s | 0.2 / 0.4 $mg$ |
| A3 forte | 6 m/s | +6 m/s | 0.3 / 0.6 $mg$ |

**B — Carico.** Massa e inerzia del corpo moltiplicate per lo stesso fattore nella scena, prima di creare il `Drone`, che quindi le legge (il controllore le conosce).

| Livello | Massa | $T/W$ | Hover a | Accelerazione angolare massima |
|---|---|---|---|---|
| B1 | ×1.5 | 2.67 | 61% di `Omega_max` | invariata |
| B2 | ×2 | 2.0 | 71% | invariata |
| B3 | ×2.5 | 1.6 | 79% | −40% |

**C — Motori lenti.** Solo $\tau_m$ in `PARAM_MOTORI`; nessun controllore lo modella. Il riferimento è la banda dell'assetto (poli a −9 rad/s).

| Livello | $\tau_m$ | Banda del motore $1/\tau_m$ | Rapporto con i poli d'assetto |
|---|---|---|---|
| nominale | 20 ms | 50 rad/s | 5.6 |
| C1 | 40 ms | 25 rad/s | 2.8 |
| C2 | 60 ms | 17 rad/s | 1.9 |
| C3 | 80 ms | 12.5 rad/s | 1.4 |

Sotto un rapporto di circa 2 la separazione tra motore e assetto non regge più: C2 e C3 sono pensati per trovare la soglia di rottura.

### 5.4 T1 — Hover

Riferimento fermo in $P$, $\psi_d = 0$. È il test delle **condizioni iniziali** e della **resistenza ai disturbi**.

| Prova | Condizione iniziale | Durata | Fasi | Cosa mette in evidenza |
|---|---|---|---|---|
| N1 | spostato di 10 cm in x | 10 s | F1 [0, 6), F4 [7, 10) | test di taratura (§2.1): le tre risposte devono quasi sovrapporsi |
| N2 | spostato di (1, 1, 0.5) m e con yaw di 90° | 10 s | F1 [0, 6), F4 [7, 10) | convergenza da errore grande: $t_c$, sovraelongazione, saturazioni |
| N3 | rilasciato a 3 m con inclinazione iniziale di 30°, 60°, 90°, 120°, 179°, **separatamente su rollio e su beccheggio** | 8 s | F1 | H3: con ZYX solo il recupero sul beccheggio passa per $\theta = \pm 90°$. Il confronto rollio/beccheggio isola l'effetto della singolarità |
| A1–A3 | ferma in $P$; vento da $t = 0$, raffica a $t_0 = 8$ s | 16 s | F2 [4, 8), F3 [8, 13), F4 [13, 16) | errore a regime sotto vento costante (nessuno ha integrale), picco e recupero dalla raffica |
| B1–B3 | come N2, drone pesante | 10 s | come N2 | convergenza con meno riserva di spinta: chi satura, chi va in windup |
| C1–C3 | come N2, motori lenti | 10 s | come N2 | oscillazioni e stabilità dell'assetto |

Cosa aspettarsi in N3 oltre a H3: il BS calcola la spinta come $\|F\|$ senza guardare l'assetto attuale. Da capovolto spinge quindi a piena potenza **verso il basso**, mentre GEO proietta la forza sull'asse del drone ($f = F \cdot R e_3$) e la riduce da solo. La quota persa dovrebbe mostrarlo.

### 5.5 T2 — Otto

Lemniscata di Gerono, $\psi_d = 0$:

$x_d = A\sin\omega t,\quad y_d = \frac{A}{2}\sin 2\omega t,\quad z_d = z_0$, con $A = 1.5$ m, $z_0 = 1.5$ m.

- Le derivate sono analitiche fino allo snap, quindi DFL e GEO ricevono un feedforward esatto. Serve una classe con la stessa interfaccia di `TraiettoriaInterpolata` (che restituisca il dict `p, v, a, j, s, psi, dpsi, ddpsi`).
- Al centro dell'otto l'accelerazione è nulla e il jerk è massimo: l'assetto passa per l'orizzontale con la massima velocità angolare. Ai lobi l'accelerazione, e quindi l'inclinazione, è massima.
- La velocità del riferimento non si annulla mai.
- 4 giri: il primo è F1, perché il drone parte fermo nel punto iniziale mentre il riferimento è già in moto; i giri 2–4 sono F2, con metriche giro per giro.

**Nominale**: scansione di $\omega$. Valori calcolati dal riferimento:

| $\omega$ [rad/s] | 0.8 | 1.0 | 1.4 | 1.8 |
|---|---|---|---|---|
| inclinazione massima | 11.7° | 18.0° | 32.5° | 46.5° |
| spinta massima | 1.02 $mg$ | 1.05 $mg$ | 1.19 $mg$ | 1.45 $mg$ |
| jerk massimo [m/s³] | 3.2 | 6.2 | 17.0 | 36.1 |
| durata di un giro | 7.9 s | 6.3 s | 4.5 s | 3.5 s |

**A, B, C** alla $\omega$ di riferimento (probabilmente 1.4). In A la raffica parte nel giro 3, quando il drone arriva al lobo $x = +A$, cioè dove la richiesta di inclinazione è già massima. F3 è il giro 3, confrontato con il giro 2.

### 5.6 T3 — Waypoint

`TraiettoriaInterpolata` con waypoint $(x, y, z, \psi)$: un quadrato di 3 m di lato, quote alternate tra 1 e 2.5 m, e lo yaw che ruota di 90° a ogni lato:

$(0, 0, 1, 0°) \to (3, 0, 2.5, 90°) \to (3, 3, 1, 180°) \to (0, 3, 2.5, 270°) \to (0, 0, 1, 360°)$

Mette insieme quello che l'otto non ha: **spigoli** (cambi di direzione bruschi), **salite e discese** (la spinta varia molto e scende nelle discese) e **yaw che ruota durante il moto**. Il drone parte fermo nel primo punto e, dopo la fine della traiettoria, resta in hover nell'ultimo per 4 s (F4).

**Nominale**: scansione di `v_media` con `min_snap`. Valori calcolati con `traiettoria.py`:

| `v_media` [m/s] | 1 | 1.5 | 2 | 2.5 |
|---|---|---|---|---|
| inclinazione massima | 8.4° | 19.9° | 37.6° | 59.6° |
| spinta minima | 0.91 $mg$ | 0.82 $mg$ | 0.71 $mg$ | 0.60 $mg$ |
| durata | 13.4 s | 8.9 s | 6.7 s | 5.4 s |

A 3 m/s si arriva a 80°: troppo, per questo la scansione si ferma a 2.5.

**Variante di regolarità** (H6). La forma resta la stessa e cambia solo quante derivate sono continue: `cubica`, `min_jerk`, `min_snap`. ⚠️ **A parità di tempi `min_snap` chiede picchi d'accelerazione più alti** (a 2 m/s: 37.6° contro 26.0° della cubica): minimizza lo snap, non il picco. Confrontare i metodi agli stessi tempi mescolerebbe due effetti. Il confronto va fatto **a parità di inclinazione massima**: `min_snap` a 2.0 m/s, `min_jerk` a 2.29 m/s, `cubica` a 2.30 m/s danno tutti 37.6°. La `lineare` è esclusa: ha gradini di velocità che nessun drone può seguire.

**A, B, C** alla `v_media` di riferimento (probabilmente 2), con `min_snap`. In A la raffica arriva durante il tratto 2, in discesa, dove la spinta è più bassa.

### 5.7 Insidie

- **GEO vicino a 180°.** L'errore d'assetto $e_R = \frac12(R_d^T R - R^T R_d)^\vee$, per una rotazione di un angolo $\vartheta$, vale $\sin\vartheta\,\hat n$. A 180° esatti è **nullo**: è un equilibrio instabile e il drone non si muoverebbe. Per questo in N3 il caso estremo è 179°. Vicino a 180° GEO parte lento, perché l'errore va come il seno: è una proprietà del metodo e va detta.
- **DFL a spinta nulla.** $J(\xi, \zeta)$ è singolare per $T = 0$. In N3, nelle discese veloci di T3 e in saturazione si può arrivare vicino. Serve una protezione (per esempio un $T$ minimo nel compensatore). La protezione cambia il metodo: va dichiarata.
- **Limiti e derivate del BS**: vedi §2.

### 5.8 Quante prove

| Test | Nominale | A + B + C | Totale per controllore |
|---|---|---|---|
| T1 | N1 + N2 + N3 (5 angoli × 2 assi) = 12 | 9 | 21 |
| T2 | 4 velocità | 9 | 13 |
| T3 | 4 velocità + 2 di regolarità = 6 | 9 | 15 |
| | | | **49** |

Con tre controllori sono circa 150 prove. Ogni prova dura 10–20 s simulati: vanno lanciate in automatico da uno script (§6), non a mano.

---

## 6. Raccolta dei risultati

### 6.1 Un log per prova

Un file per prova: `risultati/<test>/<scenario>/<controllore>.npz` (per esempio `risultati/T2/B2/GEO.npz`), con:

- **serie temporali**: $t$, $p$, $p_d$, $v$, $v_d$, $\psi$, $\psi_d$, $R$, $\omega$, $T$ e $\tau$ comandati e realizzati, $\Omega_{cmd}$, $\Omega$, flag di saturazione, $t_{calc}$, forza dell'aria; più le grandezze interne per la diagnosi ($\eta_d$ o $R_d$, $\zeta$, funzioni di Lyapunov);
- **metadati**: controllore, test, scenario e livello, guadagni, parametri di modello e del drone, `DT`, hash del commit git, data.

I log grezzi restano fuori da git (`.gitignore`); le tabelle riassuntive invece si salvano nel repo.

### 6.2 Simulazione e analisi separate

Lo script che lancia le prove salva solo i log. Un secondo script legge i log e calcola metriche, tabelle e grafici. Così si può cambiare una metrica o la soglia $\varepsilon$ senza rifare le simulazioni in CoppeliaSim.

**Verifica di determinismo** (una volta, all'inizio): lanciare due volte la stessa prova. Se i log coincidono, basta una prova per configurazione. Altrimenti servono più ripetizioni.

### 6.3 Tabelle

Una tabella per prova e fase, con i controllori in colonna; il valore migliore in grassetto, "—" se la metrica non si applica o la prova è fallita. Modello:

**T2-A2 — Otto a ω = 1.4, vento medio: giro 3 (raffica) contro giro 2**

| Metrica | BS | DFL | GEO |
|---|---|---|---|
| RMSE giro 2 / giro 3 [m] | | | |
| Picco $\|e_p\|$ [m] | | | |
| Massimo $\lvert e_x \rvert$ [m] | | | |
| $t_r$ ($\varepsilon$ = 5 cm) [s] | | | |
| RMS $\|\tau\|$ [N·m] | | | |
| Saturazione [%] | | | |

Più una **mappa di sintesi 3×4**: righe = test, colonne = N, A, B, C, e in ogni cella il controllore migliore con il distacco dal secondo. È la pagina da cui parte la discussione nella relazione.

### 6.4 Quando una differenza è rilevante

Da fissare **prima** di guardare i risultati. Una differenza si evidenzia se:

1. supera il 10% del valore migliore;
2. supera 5 mm in assoluto (sotto quel livello non ha senso fisico);
3. supera la variabilità tra prove ripetute o tra giri (su T2: intervalli media ± 2σ sui giri che non si sovrappongono).

Così non si commentano differenze che sono rumore.

### 6.5 Grafici

- **Curve di difficoltà**: RMSE in funzione dell'inclinazione massima richiesta, una curva per controllore. T2 e T3 stanno sullo stesso asse x, quindi si possono mettere nello stesso grafico.
- **Curve di degrado**, una per colonna A, B, C: la metrica principale in funzione del livello, una curva per controllore, con la soglia di rottura segnata.
- **T1-N3**: successo/fallimento e quota persa in funzione dell'inclinazione iniziale, rollio e beccheggio affiancati.
- $\|e_p\|(t)$ dei tre controllori sovrapposti, con le fasi ombreggiate e la linea di $\varepsilon$: mostra *dove* i metodi si separano.
- Traiettorie 3D sovrapposte al riferimento, per le prove più significative.
- Boxplot di $t_{calc}$.

---

## 7. Organizzazione del codice (proposta)

```
comune/
  drone.py              impianto: stato, mixer, motori
  traiettoria.py        TraiettoriaInterpolata + otto analitico + hover, stessa interfaccia
  disturbi.py           modello dell'aria (attrito, vento, raffica)
  scenari.py            test × scenario × livello: riferimento, posa iniziale, modifiche al drone, aria, durata, fasi
  metriche.py           metriche da un log
Backstepping_control/   controllore_backstepping.py
DFL_control/            controllore_dfl.py
Geometric_control/      controllore_geometrico.py
esperimenti/
  campagna.py           test × scenari × controllori → log
  analisi.py            log → tabelle (markdown/CSV) e grafici
risultati/
```

**Modifiche a `Drone`**: flag di saturazione; $T, \tau$ realizzati; motori inizializzati alla velocità di hover; per B, massa e inerzia scalate nella scena prima della costruzione.

**Interfaccia comune dei controllori.** È il prerequisito per un unico script che lanci tutte le prove:

| Elemento | Contenuto |
|---|---|
| costruttore | $m$, $J$, guadagni, limiti, `dt` |
| `reset()` | azzera lo stato interno: memorie del filtro del BS, compensatore $\zeta$ del DFL |
| `calcola(stato, rif)` | `stato` contiene $p, v, R, \eta, \omega$ e ognuno prende ciò che gli serve; `rif` è il dict della traiettoria; restituisce $(T, \tau)$ |
| grandezze interne | per il log di diagnosi: $\eta_d$ o $R_d$, Lyapunov, $\zeta$ |

Oggi il BS ha `calcola(p, v, eta, omega, rif)` e non riceve $R$: passare a un unico argomento `stato` è una modifica piccola.

---

## 8. Piano di lavoro

| Fase | Attività | Finita quando |
|---|---|---|
| 0 | Infrastruttura: interfaccia comune, modifiche a `Drone`, otto analitico, modello dell'aria, scalatura di massa e inerzia, log, metriche, verifica di determinismo. Collaudo con il BS, che esiste già | il BS gira su T1-N1 e T1-N2 da `campagna.py` e `analisi.py` produce le tabelle |
| 1 | DFL | T1-N1 e T2 a ω = 0.8 convergono |
| 2 | GEO | idem |
| 3 | Taratura a poli uguali | in T1-N1 le tre risposte sono quasi sovrapposte |
| 4 | Campagna nominale: T1-N, scansioni di T2 e T3, variante di regolarità | curve di difficoltà; scelta di $\omega$ e `v_media` di riferimento |
| 5 | Campagna degli scenari A, B, C sui tre test | mappa 3×4 e curve di degrado |
| 6 | Miglioramenti scelti da §9 | — |
| 7 | Analisi e relazione | ogni ipotesi H1–H7 è confermata o smentita con un numero o un grafico |

---

## 9. Suggerimenti di miglioramento

**Sul metodo di confronto**

1. **PD come riferimento** (slide 09, pag. 4). Costa poco ed è già nelle slide. Risponde alla domanda "quanto guadagnano i metodi non lineari rispetto a un PD?". Senza un termine di paragone, "il BS ha un RMSE di 3 cm" dice poco.
2. **Sensibilità ai guadagni**: ripetere T1-N2 e T2 alla $\omega$ di riferimento con tutti i guadagni ×0.7 e ×1.3. Se la classifica cambia, le conclusioni dipendono dalla taratura, e va scritto.
3. **Parità di sforzo** come secondo criterio di taratura: tarare i tre per lo stesso RMS di $\|\tau\|$ e confrontare gli errori.
4. **Massa non nota** come controprova di B: il drone B1 con il controllore che usa la massa nominale. Separa l'effetto "meno agile" dall'effetto "modello sbagliato" (errore a regime, nessuno ha integrale).
5. **Monte Carlo** su B e C: estrarre massa, $\tau_m$ e vento da intervalli invece dei tre livelli fissi.
6. **Un secondo motore fisico**: CoppeliaSim offre Bullet, ODE, Newton e MuJoCo. Ripetere qualche prova con un altro motore verifica che le conclusioni non dipendano dal simulatore.

**Sui controllori**

7. **Derivate analitiche dell'assetto desiderato nel BS**, da jerk e snap invece che filtrate (§2).
8. **Spinta proiettata nel BS** ($T = F \cdot R e_3$, come in GEO): dovrebbe ridurre la quota persa in T1-N3. Da confrontare come variante, non al posto della versione base.
9. **Azione integrale** (o stima del disturbo) come variante di tutti e tre. In T1-A separa la parte di errore a regime dovuta al metodo da quella dovuta alla mancanza dell'integrale.
10. **Mixer con priorità all'assetto**: quando i motori saturano, sacrificare la spinta e non le coppie (come fanno gli autopiloti reali). Cambia l'impianto per tutti allo stesso modo; interessante su B3.
11. **Compensazione della dinamica dei motori**: nessuno dei tre la modella; per il DFL è l'errore più vistoso (colonna C).

**Prove in più**

12. **Yaw tangente al percorso sull'otto** ($\psi_d = \mathrm{atan2}(\dot y_d, \dot x_d)$): aggiunge l'accoppiamento tra yaw e traslazione.
13. **Looping verticale** come dimostrazione finale: fa attraversare $\theta = 90°$, quindi BS e DFL non possono farlo per costruzione, GEO sì (è la tesi di `LeLoMc10.pdf`). A velocità costante con $v^2/r = 1.5\,g$ la spinta resta tra 0.5 e 2.5 $mg$, dentro i limiti dei motori.
14. **Stima dello stato** (lezione 06): IMU simulata + filtro complementare o di Mahony al posto dello stato ideale. Molto lavoro: solo se resta tempo.

**Limiti da dichiarare nella relazione**

- La taratura condiziona tutto: mitigata da poli uguali + sensibilità, non eliminata.
- I tempi di calcolo misurano l'implementazione Python, non l'algoritmo.
- Un solo modello di drone e un solo simulatore (salvo il punto 6).
- CoppeliaSim non modella l'aerodinamica del rotore (drag delle eliche, flappaggio, effetto suolo): c'è solo quella aggiunta con il modello dell'aria.

---

## 10. Decisioni

**Prese**

- [x] Niente replica delle slide 09: tre riferimenti nuovi (hover, otto, waypoint)
- [x] Per ogni riferimento, nominale + tre scenari che cambiano una condizione alla volta: aria, carico, motori
- [x] B varia solo la massa (con l'inerzia); `Omega_max` resta fisso in tutte le prove

**Aperte**

- [ ] Livelli di A, B, C (tabelle in §5.3) e coefficiente $c_a$
- [ ] Geometria dei waypoint di T3
- [ ] Limiti comuni: nessun limite d'inclinazione, o lo stesso per tutti (60°?)
- [ ] Derivate dell'assetto desiderato ($\dot\phi_d, \ddot\phi_d, \dot\theta_d, \ddot\theta_d$ nel BS; $\Omega_d, \dot\Omega_d$ in GEO): filtrate o analitiche, con la stessa scelta per entrambi
- [ ] Poli di riferimento: quelli dei guadagni attuali del BS o altri?
- [ ] Soglia $\varepsilon$ (5 cm?)
- [ ] Includere il PD come riferimento?
- [ ] Velocità di riferimento per gli scenari di T2 e T3: si decide dopo la campagna nominale (fase 4)
