# Confronto BS / DFL / GEO — piano passo passo

Guida **senza codice** per realizzare il progetto, dalla taratura alla presentazione: il codice lo scrivi tu.
I numeri attesi sono nascosti nei blocchi **«Verifica»**: prova a ricavarli prima di aprirli.

**Documenti collegati**

- `doc/definizione_proposta/definizione_proposta.tex` — proposta approvata: è il riferimento vincolante
- `planning.md` — dettagli operativi: livelli degli scenari, metriche, fasi del volo, ipotesi H1–H7
- `doc/Taratura del backstepping.html` — progetto dei poli del BS sulla carta di Nichols, con grafici interattivi
- `Backstepping_control/STEPS_backstepping.md` — derivazione del backstepping

**Il filo del piano.** La FASE 0 spiega da dove vengono tutti i numeri della taratura. Le fasi 1–4 costruiscono il codice e verificano che la taratura sia davvero la stessa per i tre controllori. Le fasi 5–8 sono le campagne di prove e l'analisi.

| Fase | Contenuto | Finita quando |
|---|---|---|
| 0 | Analisi lineare in hover: da dove vengono poli e guadagni | `analisi/lineare.py` riproduce i numeri dei blocchi Verifica |
| 1 | La taratura nel codice: i poli come unica fonte | i guadagni dei tre controllori escono da una sola funzione |
| 2 | Organizzazione del codice per la combinatoria delle prove | il BS gira su T1-N1 e T1-N2 lanciato da un main, con un log per prova |
| 3 | DFL e GEO | T1-N1 e T2 a ω = 0.8 convergono |
| 4 | Verifica della taratura | T1-N1 coincide con la previsione lineare per tutti e tre |
| 5 | Campagna nominale | curve di difficoltà e difficoltà di riferimento scelta |
| 6 | Campagna degli scenari A, B, C | mappa 3×4 e curve di degrado |
| 7 | *(facoltativa)* Sensibilità alla taratura | classifica a κ = 0.7, 1, 1.3 e decisione per la presentazione |
| 8 | Analisi e presentazione | ogni ipotesi H1–H7 ha un verdetto |

---

## FASE 0 — Da dove vengono i numeri

L'idea in una riga: **i guadagni si scelgono in modo che i tre controllori, linearizzati attorno all'hover, abbiano gli stessi poli in anello chiuso.** In hover i tre diventano così la stessa legge lineare, e ogni differenza nelle prove viene dalla costruzione dei metodi.

Tutto ciò che segue si fa su carta e con uno script offline, senza CoppeliaSim.

### 0.1 Il modello in hover

**Esercizio 1** — Linearizza il modello attorno all'hover: $p = p_0$, $v = 0$, $\phi = \theta = \psi = 0$, $\omega = 0$, $T = mg + \delta T$. Scrivi le equazioni di $x, y, z, \phi, \theta, \psi$.

<details>
<summary>Verifica</summary>

$$\ddot x = g\,\theta,\qquad \ddot y = -g\,\phi,\qquad \ddot z = \frac{\delta T}{m}$$

$$\ddot\phi = \frac{\tau_\phi}{I_x},\qquad \ddot\theta = \frac{\tau_\theta}{I_y},\qquad \ddot\psi = \frac{\tau_\psi}{I_z}$$

Ogni asse laterale è una **catena di 4 integratori**: $x \to \dot x \to \theta \to \dot\theta$, con un solo ingresso ($\tau_\theta$). Quota e yaw sono doppi integratori. Gli assi sono disaccoppiati.
</details>

**Esercizio 2** — Aggiungi le due dinamiche che nessun controllore conosce: i motori e il campionamento.

<details>
<summary>Verifica</summary>

- **Motori**: ogni ingresso ($T$ e le tre coppie) passa per $\dfrac{1}{\tau_m s + 1}$. Il mixer è lineare in $\Omega^2$, e linearizzando $\Omega^2$ attorno a $\Omega_h$ il primo ordine sulla velocità diventa un primo ordine su $T$ e $\tau$.
- **Campionamento**: il comando resta costante per $T_c = 5$ ms, che equivale circa a un ritardo di mezzo passo, $e^{-sT_c/2}$.

Massa e inerzia si semplificano, perché i tre controllori le cancellano. L'analisi dipende solo da $g$, dai poli, da $\tau_m$ e da $T_c$.
</details>

### 0.2 Le tre leggi linearizzate

Notazione: per ogni canale il polinomio del secondo ordine è $s^2 + \beta s + \alpha$, con $\beta = 2\zeta\omega_n$ e $\alpha = \omega_n^2$. Sull'asse x: posizione $(\beta_p, \alpha_p)$, assetto $(\beta_a, \alpha_a)$. Nel documento HTML: $a_1 = \beta_p$, $a_0 = \alpha_p$, $b_1 = \beta_a$, $b_0 = \alpha_a$.

**Esercizio 3** — **BS** sull'asse x. Linearizza l'inversione ($\theta_d \approx \mu_x / g$) e il controllo d'assetto. Ricava il polinomio caratteristico in due casi: derivate $\dot\theta_d, \ddot\theta_d$ **nulle** (`tau_filtro=None`) e **esatte**.

<details>
<summary>Verifica</summary>

Il singolo canale del BS dà $s^2 + (k_1 + k_2)s + (1 + k_1 k_2)$.

Derivate nulle:
$$s^4 + \beta_a s^3 + \alpha_a s^2 + \alpha_a\beta_p\, s + \alpha_a\alpha_p$$

Derivate esatte: l'errore d'assetto $e_\theta = \theta_d - \theta$ evolve da solo ($\ddot e_\theta + \beta_a \dot e_\theta + \alpha_a e_\theta = 0$) e la posizione ne è solo forzata. Il sistema è triangolare a blocchi, quindi
$$(s^2 + \beta_p s + \alpha_p)(s^2 + \beta_a s + \alpha_a)$$

**Solo con le derivate esatte** i poli sono quelli scelti per i due canali.
</details>

**Esercizio 4** — **DFL** sull'asse x. La legge impone $x^{(4)} = x_d^{(4)} + c_3 e^{(3)} + c_2 \ddot e + c_1 \dot e + c_0 e$. In hover $\ddot x = g\theta$ e $x^{(3)} = g\dot\theta$. Qual è il polinomio?

<details>
<summary>Verifica</summary>

$$s^4 + c_3 s^3 + c_2 s^2 + c_1 s + c_0$$

Nel modello nominale vale **ovunque**, non solo in hover: la cancellazione è esatta (salvo $T = 0$). Per avere gli stessi poli del BS, i $c_i$ sono i coefficienti del prodotto dei due polinomi.
</details>

**Esercizio 5** — **GEO**. In hover $e_R \approx \eta - \eta_d$ e $f \approx m(g + \ddot z)$. Ricava i polinomi di posizione e assetto.

<details>
<summary>Verifica</summary>

$$s^2 + \frac{k_v}{m}s + \frac{k_x}{m}\qquad\qquad s^2 + \frac{k_\Omega}{J_i}s + \frac{k_R}{J_i}$$

È la stessa struttura a cascata del BS. Con $\Omega_d, \dot\Omega_d$ esatte, il polinomio dell'asse x è il prodotto dei due; con $\Omega_d = 0$ è quello "a derivate nulle" del BS.
</details>

### 0.3 Perché "stessi poli" vuol dire "stessa legge"

**Esercizio 6** — Mostra che, in hover, ognuna delle tre leggi sull'asse x diventa $\tau_\theta = -K\,(x, \dot x, \theta, \dot\theta)^T$. Perché due leggi con lo stesso polinomio caratteristico hanno lo stesso $K$?

<details>
<summary>Verifica</summary>

La catena ha un solo ingresso ed è controllabile. Per un sistema così il polinomio caratteristico fissa $K$ in modo unico (formula di Ackermann). Se i polinomi coincidono, coincidono le leggi.

Conseguenze:

- su **x, y, ψ** i tre controllori linearizzati sono **la stessa legge**, non tre leggi simili;
- su **z** no: il DFL ha due stati in più, il compensatore $\zeta = (T, \dot T)$ (vedi 0.7);
- l'ipotesi H1 vale per costruzione in hover. T1-N1 quindi **verifica l'implementazione**, non confronta i metodi;
- nelle altre prove le differenze vengono solo dai termini non lineari e dall'asse z, cioè dalla costruzione.
</details>

> ⚠️ Vale solo se le derivate dell'assetto desiderato sono esatte o ben approssimate (vedi 0.6). "Derivate analitiche da jerk e snap" **non basta**: in hover jerk e snap del riferimento sono nulli, mentre $R_d$ cambia per il termine di retroazione $-k_x e_x - k_v e_v$. Per derivarlo servono accelerazione e jerk *del drone*, dal modello: $a = \frac{T}{m}Re_3 - g e_3$, $j = \frac{1}{m}(\dot T R e_3 + T R\hat\omega e_3)$.

### 0.4 Quanto veloci: le specifiche sulla carta di Nichols

L'assegnamento dei poli dice **dove** mettere i poli, ma lavora sul modello che il controllore conosce, dove si potrebbero rendere veloci quanto si vuole. Il limite viene da motori e campionamento, che tolgono fase, e ne tolgono di più alle frequenze alte. La carta di Nichols disegna $L(j\omega)$ (modulo in dB contro fase) e da lì si leggono margine di fase, margine di guadagno e picco di risonanza. Il ritardo vi compare come un semplice sfasamento $-\omega T_c/2$. È il passo in cui le specifiche di robustezza entrano nella taratura. I dettagli e i grafici sono nel documento HTML.

**Specifica scelta**: margine di fase ≥ 50°, con motore e campionamento inclusi.

**Esercizio 7** — Vincolo del BS: per quali $(\omega_n, \zeta)$ esistono $k_1, k_2$ reali e positivi?

<details>
<summary>Verifica</summary>

$k_1$ e $k_2$ sono le radici di $k^2 - \beta k + (\alpha - 1) = 0$. Sono reali e positive se e solo se

$$\omega_n > 1\ \text{rad/s},\qquad \zeta \ge \sqrt{1 - 1/\omega_n^2}$$

Quindi i poli complessi del BS hanno sempre **parte immaginaria ≤ 1 rad/s**. DFL e GEO non hanno questo vincolo: i poli comuni vanno scelti dentro la zona ammessa dal BS.

- Assetto con $\omega_n = 9$: serve $\zeta \ge 0.994$, quindi $\zeta = 1$ (polo doppio).
- Posizione con $\omega_n = 1.58$: serve $\zeta \ge 0.775$, e $\zeta = 0.79$ va bene.
</details>

**Esercizio 8** — Anello d'assetto: $L(s) = \dfrac{(\beta_a s + \alpha_a)\,e^{-sT_c/2}}{s^2(\tau_m s + 1)}$ con $\zeta = 1$. Fino a quale $\omega_n$ il margine di fase resta sopra 50°?

<details>
<summary>Verifica</summary>

Fino a $\omega_n \approx 10.8$ rad/s. Con $\omega_n = 9$: margine di fase 54° e margine di guadagno 26 dB. Quindi $\beta_a = 18$, $\alpha_a = 81$ e $k = (8, 10)$, i guadagni di `main_hover.py`.
</details>

**Esercizio 9** — Anello di posizione: aprilo dove entra $\theta_d$ e poi all'attuatore, dove sta il motore. Confronta derivate nulle e derivata filtrata ($\tau_f = 20$ ms).

<details>
<summary>Verifica</summary>

| | Derivate nulle | Derivata filtrata, $\tau_f$ = 20 ms |
|---|---|---|
| margine di fase aprendo su $\theta_d$ | 37° | ≈ 70° |
| margine di fase all'attuatore | — | 43° |

Il margine che conta è il più piccolo: 43° con la derivata filtrata.
</details>

**I poli nominali** che escono da 0.4:

| Canale | Poli | $\omega_n$, $\zeta$ | BS $(k_1, k_2)$ |
|---|---|---|---|
| x, y | $-1.25 \pm 0.97j$ | 1.58, 0.79 | (1, 1.5) |
| z | $-2.5 \pm 0.87j$ | 2.65, 0.945 | (2, 3) |
| φ, θ | −9 doppio | 9, 1 | (8, 10) |
| ψ | −5 doppio | 5, 1 | (4, 6) |
| z, coppia in più del DFL | **da decidere** (0.7) | | — |

### 0.5 Dai poli ai guadagni dei tre controllori

**Esercizio 10** — Calcola $(\beta, \alpha)$ di ogni canale e da lì i guadagni del BS, i coefficienti del DFL e le matrici di GEO.

<details>
<summary>Verifica</summary>

- **BS**: $k_{1,2} = \dfrac{\beta \pm \sqrt{\beta^2 - 4(\alpha - 1)}}{2}$. Con un polo doppio in $-\lambda$: $k = \lambda \pm 1$.
- **DFL**:
  - x, y: coefficienti di $(s^2 + 2.5s + 2.5)(s^2 + 18s + 81)$, cioè $c_3 \dots c_0 = 20.5,\ 128.5,\ 247.5,\ 202.5$;
  - ψ: $(10, 25)$;
  - z: coefficienti di $(s^2 + 5s + 7)$ per la coppia in più.
- **GEO**, con matrici diagonali:
  - $K_x = m\,\mathrm{diag}(2.5,\ 2.5,\ 7)$, $K_v = m\,\mathrm{diag}(2.5,\ 2.5,\ 5)$;
  - $K_R = \mathrm{diag}(81 I_x,\ 81 I_y,\ 25 I_z)$, $K_\Omega = \mathrm{diag}(18 I_x,\ 18 I_y,\ 10 I_z)$.
</details>

> ⚠️ **L'invariante sono i poli, non i guadagni.** Nello scenario B massa e inerzia cambiano e il controllore le conosce. BS e DFL mantengono i poli da soli, perché moltiplicano per $m$ e $J$. GEO li mantiene solo se $K$ viene ricalcolato con la massa nuova.

> ⚠️ Con guadagni **scalari**, come nell'articolo di Lee, GEO non può dare a z poli diversi da x e y, né allo yaw poli diversi da rollio e beccheggio. Servono matrici diagonali. La dimostrazione di Lee è scritta per guadagni scalari: la scelta va dichiarata.

### 0.6 Le derivate dell'assetto desiderato

**Esercizio 11** — Gradino di 10 cm su x, con il motore. Confronta $x(t)$ del BS con le derivate esatte, filtrate e nulle.

<details>
<summary>Verifica</summary>

| Derivate di $\theta_d$ (BS) o $\Omega_d$ (GEO) | Poli effettivi (senza motore) | Distanza massima dalla risposta esatta |
|---|---|---|
| esatte | −9, −9, −1.25 ± 0.97j | 0 |
| filtrate, $\tau_f$ = 20 ms | praticamente quelli | 0.4 mm |
| nulle | −12.8, −1.9, −1.64 ± 2.37j | **19 mm su 100** |

La separazione tra le scale temporali non basta: anche con rapporto 6 tra assetto e posizione, con le derivate nulle si sbaglia del 19%.
</details>

**Decisione**: la stessa scelta per BS e GEO. Proposta: derivata filtrata con $\tau_f$ = 20 ms (almeno 4 passi di controllo), dichiarata nella relazione.

### 0.7 L'asse z del DFL

**Esercizio 12** — Nel DFL z è del quarto ordine: $(s^2 + 5s + 7)$ per una coppia di poli in più da scegliere. Confronta un gradino di 0.5 m su z con quello di BS e GEO.

<details>
<summary>Verifica</summary>

Con il motore, BS e GEO arrivano entro il 5% in 1.64 s.

| Poli in più del DFL | Entro il 5% | Distanza massima da BS/GEO |
|---|---|---|
| −9 doppio | 1.88 s | 109 mm (22%) |
| −15 doppio | 1.78 s | 70 mm |
| −25 doppio | 1.72 s | 43 mm (9%) |
</details>

**Decisione**: −9 doppio (come nel planning) oppure −20/−25 doppio (più vicino a BS/GEO; controllare in N2 che la spinta non saturi). In ogni caso z è l'asse in cui "stessi poli" non vuol dire "stessa legge": va dichiarato.

### 0.8 Cosa prevede l'analisi per lo scenario C

**Esercizio 13** — Calcola il margine di fase all'attuatore dell'asse x al variare di $\tau_m$, con campionamento e derivata filtrata.

<details>
<summary>Verifica</summary>

| $\tau_m$ | 20 ms (N) | 40 ms (C1) | 60 ms (C2) | 80 ms (C3) |
|---|---|---|---|---|
| margine di fase all'attuatore | 43° | 28° | 18° | 11° |
| margine di fase del solo anello d'assetto | 54° | 40° | 30° | 23° |

L'anello lineare diventa instabile verso $\tau_m \approx 120$ ms.
</details>

Due conseguenze da scrivere **prima** di simulare:

- su x, y e ψ la legge lineare è la stessa per tutti e tre, quindi in hover i margini perdono lo stesso valore per tutti. Differenze in C possono venire solo da z e dalle non linearità (H5);
- C3 = 80 ms probabilmente non trova la soglia di rottura in hover. **Decisione**: lasciare i livelli o aggiungerne uno verso 100–120 ms.

### 0.9 Lo script `analisi/lineare.py`

Riproduce tutta la FASE 0, senza CoppeliaSim.

| | Contenuto |
|---|---|
| Ingressi | i poli (dalla FASE 1), $\tau_m$, $T_c$, $\tau_f$ |
| Uscite | polinomi; autovalori dell'anello chiuso con e senza motore; margini e carta di Nichols; risposte al gradino previste (x a 10 cm, z a 0.5 m) |
| Strumenti | modello in spazio di stato; autovalori con `numpy.linalg.eigvals`; risposta con l'esponenziale di matrice; ritardo come $e^{-j\omega T_c/2}$ in frequenza |

**Finita quando** riproduce i numeri dei blocchi Verifica. La carta di Nichols e le previsioni per C sono candidate per la presentazione.

---

## FASE 1 — La taratura nel codice

Principio: **una sola fonte di verità, i poli.** I guadagni si calcolano, non si scrivono mai a mano. Così è impossibile che un controllore resti con guadagni vecchi.

`comune/taratura.py`:

| Elemento | Ingressi | Uscita |
|---|---|---|
| `POLI` | — | $(\omega_n, \zeta)$ per canale (x/y, z, φ/θ, ψ), coppia in più del DFL su z, $\tau_f$ |
| `polinomi(poli)` | poli | $(\beta, \alpha)$ per canale |
| `guadagni_bs(poli)` | poli | dict nel formato `GUADAGNI` di `main_hover.py`; errore se un canale è fuori dal vincolo dell'Esercizio 7 |
| `guadagni_dfl(poli)` | poli | $c_3 \dots c_0$ per x, y e z, $(c_1, c_0)$ per ψ |
| `guadagni_geo(poli, m, J)` | poli, massa e inerzia lette dalla scena | $K_x, K_v, K_R, K_\Omega$ diagonali |
| `scala_poli(poli, kappa)` | poli, $\kappa$ | poli con $\omega_n \times \kappa$ (serve solo nella FASE 7; $\kappa = 1$ di default) |

**Chi la chiama**: lo svolgimento di ogni prova (FASE 2), *dopo* aver creato il `Drone`, perché GEO ha bisogno di $m$ e $J$ della scena.

**Test**: `guadagni_bs(POLI)` deve restituire esattamente i `GUADAGNI` di oggi, cioè (1, 1.5), (2, 3), (8, 10), (4, 6).

---

## FASE 2 — Organizzazione del codice per la combinatoria delle prove

Qui poli e guadagni sono **fissati** (κ = 1). Il problema è lanciare circa 150 prove in modo automatico, ripetibile e confrontabile.

### 2.1 Da dove si parte

La proposta approvata fissa: classe del drone, classi dei controllori, classe delle traiettorie, e **un main per ogni coppia riferimento–scenario**, cioè 12 main. Ogni main esegue tutti i livelli per tutti e tre i controllori. I risultati vanno in `output/`, una sottocartella per prova, con nome dato dalla quaterna (per esempio `output/T2_A2_DFL`).

> ⚠️ Il §6.1 del planning usa `risultati/<test>/<scenario>/<controllore>.npz`. Vale la proposta: `output/<id>/`.

### 2.2 Struttura

```
comune/
  drone.py           impianto (modifiche in 2.5)
  traiettoria.py     + otto analitico e hover, stessa interfaccia di TraiettoriaInterpolata
  aria.py            attrito, vento costante, raffica 1-coseno
  taratura.py        poli -> guadagni (FASE 1)
  scenari.py         catalogo delle prove: (riferimento, scenario, livello) -> descrizione
  esecuzione.py      esegue UNA prova: (descrizione, controllore) -> log
  registro.py        salva e rilegge un log
Backstepping_control/controllore_backstepping.py
DFL_control/controllore_dfl.py
Geometric_control/controllore_geometrico.py
prove/
  T1_N.py  T1_A.py  T1_B.py  T1_C.py
  T2_N.py  T2_A.py  T2_B.py  T2_C.py
  T3_N.py  T3_A.py  T3_B.py  T3_C.py
analisi/
  lineare.py         FASE 0, senza CoppeliaSim
  metriche.py        log -> metriche per fase
  tabelle.py         metriche -> tabelle e grafici
output/              una sottocartella per prova, fuori da git
```

### 2.3 Il principio: separare il *che cosa* dal *come*

- **Che cosa**: la descrizione di una prova è un **dato**, non codice. Sta in `scenari.py`.
- **Come**: una sola funzione in `esecuzione.py` esegue qualunque prova. I 12 main restano corti: elencano livelli e controllori e chiamano quella funzione.
- **Chi non sa cosa**: i controllori non sanno in che prova sono; il drone non sa quale controllore lo comanda; la descrizione della prova non contiene guadagni.

Il vantaggio: il ciclo di controllo è scritto una volta sola. Una correzione vale per tutte le prove, e nessuna differenza tra controllori può nascere da due cicli scritti in modo diverso.

### 2.4 La descrizione di una prova

| Campo | Esempio | Chi lo usa |
|---|---|---|
| `id` | `T2_A2` | nome della cartella (più il controllore) |
| riferimento | otto con ω = 1.4 | ciclo, a ogni passo |
| posa iniziale | posizione, rotazione (le inclinazioni di N3) | ripristino della scena |
| modifiche al drone | fattore di massa e inerzia (B), $\tau_m$ (C) | prima di creare il `Drone` |
| aria | $c_a$, vento, raffica ($t_0$, $T_g$, $\Delta v$) | ciclo, a ogni passo |
| durata | 10 s | ciclo |
| fasi del volo | F1 = [0, 6), F4 = [7, 10) | solo l'analisi |
| criteri di fallimento | errore > 2 m, inclinazione > 80°, contatto col suolo, NaN | ciclo |

`scenari.py` offre una funzione che, data la coppia (riferimento, scenario), restituisce l'elenco delle descrizioni, una per livello. I valori sono quelli del planning (§5.3–§5.6).

### 2.5 Lo svolgimento di una prova

1. **Scena pulita**: simulazione ferma, scena alle condizioni iniziali.
2. **Modifiche al drone**: massa e inerzia scalate nella scena (B); $\tau_m$ nei `PARAM_MOTORI` (C).
3. **Drone**: si crea, e legge $m$ e $J$ dalla scena.
4. **Guadagni**: `taratura.py` con $m$ e $J$ appena letti.
5. **Controllore**: si crea e si chiama `reset()`.
6. **Condizioni iniziali**: posa della prova, drone fermo, motori già a $\Omega_h$; per il DFL $\zeta(0) = (mg, 0)$.
7. **Ciclo**: stato → riferimento → *cronometro* → `calcola` → *cronometro* → aria → comando ai motori → log → `sim.step()`. A ogni passo si controllano i criteri di fallimento.
8. **Fine**: simulazione ferma, log e metadati salvati, massa e inerzia della scena **rimesse ai valori originali**.

Modifiche necessarie a `Drone` (planning §7): flag di saturazione; $T$ e $\tau$ realizzati, oltre a quelli comandati; motori inizializzati a $\Omega_h$; $\tau_m$ modificabile.

> ⚠️ Se nel passo 8 la massa non viene ripristinata, la prova successiva parte con il drone pesante: un errore silenzioso che contamina tutto. Conviene salvare $m$ e $J$ originali all'inizio di ogni prova e verificarli alla fine.

> ⚠️ **Determinismo**: lancia due volte la stessa prova (T1-N1 con il BS). Se i log coincidono, basta una prova per configurazione; altrimenti servono ripetizioni (planning §6.2).

### 2.6 La combinatoria

Ogni main fissa riferimento e scenario. Ripete poi due cicli: livelli all'esterno, controllori all'interno. Se la campagna si interrompe, i livelli completati hanno così tutti e tre i controllori, e restano confrontabili.

**Ripresa**: se `output/<id>_<controllore>/` esiste già con lo stesso commit e gli stessi parametri, la prova si salta. Un'opzione permette di rifarla.

**Nomi dei livelli** (l'`id` completo è `<riferimento>_<scenario><livello>_<controllore>`):

| Main | Livelli | Esempio |
|---|---|---|
| `T1_N` | N1 (10 cm), N2 (1, 1, 0.5 m e 90° di yaw), N3-R030 … N3-R179 (rollio), N3-P030 … N3-P179 (beccheggio) | `T1_N3-P120_GEO` |
| `T2_N` | N1 … N4 = ω 0.8, 1.0, 1.4, 1.8 | `T2_N3_BS` |
| `T3_N` | N1 … N4 = `v_media` 1, 1.5, 2, 2.5 (`min_snap`); N5 = `min_jerk` a 2.29; N6 = `cubica` a 2.30 | `T3_N5_DFL` |
| `T1_A/B/C` | 1, 2, 3 (dal gradino grande, come N2) | `T1_B2_DFL` |
| `T2_A/B/C`, `T3_A/B/C` | 1, 2, 3 alla difficoltà di riferimento | `T2_A2_DFL` |

La difficoltà di riferimento di T2 e T3 si sceglie nella FASE 5. In `scenari.py` è un parametro, non un numero scritto in più punti.

**Conteggio**:

| Test | Nominale | A + B + C | Per controllore |
|---|---|---|---|
| T1 | 12 | 9 | 21 |
| T2 | 4 | 9 | 13 |
| T3 | 6 | 9 | 15 |
| | | | **49, quindi 147 in tutto** |

### 2.7 Il log di una prova

In `output/<id>/`:

- **serie temporali**: quelle del planning §6.1 ($t$, $p$, $p_d$, $R$, $\omega$, $T$ e $\tau$ comandati e realizzati, $\Omega$, saturazione, $t_{calc}$, aria, grandezze interne);
- **metadati**: quaterna, **poli e κ**, guadagni effettivi del controllore, $m$, $J$, `PARAM_MOTORI`, passo, parametri dell'aria, commit git, data, esito.

Si salvano sia i poli sia i guadagni: i poli dicono che cosa si voleva, i guadagni che cosa è stato usato davvero. Se un giorno non corrispondono, il log lo mostra.

### 2.8 Collaudo con il BS

Prima di scrivere DFL e GEO: `T1_N` con il solo BS, livelli N1 e N2. **Finita quando** i log sono in `output/`, la ripresa funziona e `analisi/metriche.py` produce la tabella di T1-N1.

---

## FASE 3 — DFL e GEO

Stessa interfaccia del BS (planning §7) e guadagni presi da `taratura.py`. Decisioni da applicare:

- derivate dell'assetto desiderato trattate come nel BS (0.6);
- coppia di poli in più del DFL su z (0.7);
- protezione del DFL a spinta nulla, da dichiarare (planning §5.7);
- limiti del BS (`angolo_max`, 90% di `T_max`): toglierli, o applicare gli stessi limiti a tutti.

**Finita quando** T1-N1 e T2 a ω = 0.8 convergono.

---

## FASE 4 — Verifica della taratura

### 4.1 Offline, prima di CoppeliaSim

Un piccolo impianto lineare in Python (le equazioni di 0.1 con il motore), comandato dalle **vere classi** dei controllori. Gradino di 10 cm su x e di 0.5 m su z. La risposta deve coincidere con quella prevista da `analisi/lineare.py`.

Trova errori di segno, di terna e di traduzione dei guadagni senza dover lanciare il simulatore.

### 4.2 T1-N1 in CoppeliaSim

Nel gradino di 10 cm l'inclinazione resta sotto 0.8°: il drone è in pieno regime lineare. Le tre risposte $x(t)$ devono coincidere tra loro e con la previsione lineare, entro una soglia fissata **prima** (proposta: 2 mm, il 2% del gradino).

Se non succede, l'errore è nell'implementazione, non nel metodo, e va corretto prima di qualsiasi confronto.

**Per la presentazione**: il grafico delle tre risposte sovrapposte alla previsione lineare dimostra che la taratura è la stessa.

---

## FASE 5 — Campagna nominale

`T1_N`, `T2_N`, `T3_N` per i tre controllori, e curve di difficoltà (errore in funzione dell'inclinazione richiesta).

**Finita quando** è scelta la difficoltà di riferimento per gli scenari: la più alta a cui tutti e tre riescono con margine. Si scrive in `scenari.py` e si fa un commit.

---

## FASE 6 — Campagna degli scenari

I 9 main degli scenari. Prima di lanciarli, metti per iscritto le previsioni della FASE 0, per esempio i margini nello scenario C, accanto alle ipotesi H1–H7.

**Finita quando** sono pronte la mappa 3×4 e le curve di degrado.

---

## FASE 7 — *(facoltativa)* Sensibilità alla taratura

**A cosa serve**: mostrare che la classifica non dipende dai poli scelti. I poli si scalano di κ = 0.7, 1, 1.3, cioè $\omega_n \times \kappa$ su tutti i canali, **per tutti e tre i controllori, DFL compreso**. Il DFL va scalato anche lui per due motivi. Se resta fermo, la parità salta. E anche per lui i poli cambiano sforzo, saturazione e margini, perché motori e saturazione non vengono cancellati.

**Esercizio 14** — Calcola le tre tarature e i loro margini con `analisi/lineare.py`.

<details>
<summary>Verifica</summary>

Guadagni del BS ($\omega_n$, $\zeta$ → $k_1, k_2$):

| Canale | κ = 0.7 | κ = 1 | κ = 1.3 |
|---|---|---|---|
| x, y | 1.11, 0.79 → (0.14, 1.61) | 1.58, 0.79 → (1, 1.5) | 2.06, **0.874** → (1.80, 1.80) |
| z | 1.85, 0.945 → (0.96, 2.55) | 2.65, 0.945 → (2, 3) | 3.44, **0.957** → (3.29, 3.29) |
| φ, θ | −6.3 doppio → (5.3, 7.3) | −9 doppio → (8, 10) | −11.7 doppio → (10.7, 12.7) |
| ψ | −3.5 doppio → (2.5, 4.5) | −5 doppio → (4, 6) | −6.5 doppio → (5.5, 7.5) |

$(\beta, \alpha)$ per DFL e GEO:

| Canale | κ = 0.7 | κ = 1 | κ = 1.3 |
|---|---|---|---|
| x, y | 1.75, 1.225 | 2.5, 2.5 | 3.59, 4.22 |
| z | 3.5, 3.43 | 5, 7 | 6.58, 11.83 |
| φ, θ | 12.6, 39.7 | 18, 81 | 23.4, 136.9 |
| ψ | 7, 12.25 | 10, 25 | 13, 42.25 |
| DFL x, y: $c_3 \dots c_0$ | 14.35, 62.96, 84.89, 48.61 | 20.5, 128.5, 247.5, 202.5 | 26.99, 225.2, 590.5, 578.3 |

Margini e gradino di 10 cm (motore, campionamento, derivata filtrata):

| | κ = 0.7 | κ = 1 | κ = 1.3 |
|---|---|---|---|
| margine di fase, anello d'assetto | 60° | 54° | 48° |
| margine di fase all'attuatore | 51° | 43° | 36° |
| margine di fase, anello di posizione | 70° | 70° | 73° |
| tempo per entrare entro 5 mm | 3.38 s | 2.37 s | 2.07 s |
| inclinazione massima | 0.4° | 0.8° | 1.4° |
</details>

> ⚠️ **A κ = 1.3 la scalatura non è pura.** Con $\omega_n = 2.06$ il BS non può tenere $\zeta = 0.79$ (Esercizio 7): $\zeta$ va alzato al minimo ammesso, 0.874 su x e y e 0.957 su z, dove $k_1 = k_2$. Gli stessi poli corretti vanno dati anche a DFL e GEO.

> ⚠️ **κ = 1.3 è sotto la specifica dei 50°.** Va bene come livello aggressivo, ma solo in nominale: insieme allo scenario C diventerebbe una prova di stabilità, non di sensibilità.

**Prove**: T1-N2 e T2 alla ω di riferimento, per 3 κ e 3 controllori. Sono 18 prove, di cui 12 nuove (κ = 1 esiste già). Nel codice κ è un parametro di `taratura.py` salvato nei metadati, e l'`id` prende un suffisso (per esempio `T2_N3_BS_k1.3`).

**Decisione per la presentazione**:

- **classifica invariata**: una tabella e una frase nella relazione; in presentazione una slide facoltativa, utile se viene chiesto quanto pesa la taratura;
- **classifica cambiata**: va in presentazione, perché limita le conclusioni: in quella prova il risultato dipende anche dalla taratura.

---

## FASE 8 — Analisi e presentazione

Metriche per fase, tabelle, mappa 3×4, curve di difficoltà e di degrado, un verdetto per ciascuna ipotesi H1–H7 (planning §4–§6).

**Candidati per la presentazione dalla FASE 0**:

1. l'argomento "stessi poli = stessa legge in hover" (0.3), che giustifica il criterio di taratura;
2. la carta di Nichols dell'assetto con la specifica dei 50° (0.4);
3. T1-N1: tre risposte sovrapposte alla previsione lineare (4.2);
4. le previsioni per lo scenario C confrontate con quanto osservato (0.8).

---

## Decisioni aperte

| Decisione | Opzioni | Dove entra |
|---|---|---|
| Derivate dell'assetto desiderato | filtrate con $\tau_f$ = 20 ms (proposta) / analitiche, includendo la parte di retroazione | 0.6, FASE 3 |
| Coppia di poli in più del DFL su z | −9 doppio (planning) / −20…−25 doppio | 0.7 |
| Livelli dello scenario C | 40, 60, 80 ms / aggiungere 100–120 ms | 0.8, `scenari.py` |
| Limiti del BS | toglierli / gli stessi per tutti | FASE 3 |
| Soglia di T1-N1 | 2 mm? | 4.2 |
| Soglia ε, PD come riferimento | planning §10 | FASE 5 |
