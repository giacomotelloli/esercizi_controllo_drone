"""
Video di una prova (opzione --video dei simula.py), salvato con il log in out/<id>_<controllore>/video.mp4.

    Ripresa        telecamera fissa aggiunta alla scena, riferimento disegnato, scia del drone, scritte; mp4 con ffmpeg
    serve_ffmpeg   controlla che ffmpeg ci sia, prima di lanciare le prove

Un fotogramma ogni PASSI passi di simulazione, cioè ogni 1/FPS s di tempo SIMULATO: il video va a velocità reale
anche se la simulazione è più lenta. Le immagini arrivano da sim.getVisionSensorImg e ffmpeg le scrive in H.264.

Riprendere non cambia la prova: telecamera e disegni non hanno massa né collisioni, e tutto si fa tra un sim.step()
e il successivo, a simulazione ferma (modo stepping). Cresce solo il tempo reale della prova, non t_calc.
La telecamera sparisce alla prova dopo, perché la scena si ricarica dal file.
"""
import shutil
import subprocess
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from comune import scenari

NOME_FILE = 'video.mp4'

LARGHEZZA, ALTEZZA = 1280, 720  # [pixel]
FPS = 25                        # fotogrammi per secondo di tempo simulato
ANGOLO_VISTA = np.radians(60)   # apertura della telecamera sul lato lungo dell'immagine [rad]
AZIMUT = np.radians(-60)        # da dove guarda la telecamera, attorno all'asse z del mondo [rad]
ELEVAZIONE = np.radians(35)     # sopra l'orizzonte [rad]
MARGINE = 0.6                   # spazio attorno alla zona del volo [m]: il drone resta dentro anche sul bordo
DISTANZA_MIN = 4.0              # [m]: nelle prove di hover la zona del volo è quasi un punto
ATTESA_FINALE = 2.0             # [s] ultimo fotogramma tenuto fermo, con l'esito della prova

COLORE_RIF = [0.15, 0.45, 0.95]     # riferimento e punto desiderato p_d
COLORE_SCIA = [0.95, 0.45, 0.10]    # traiettoria del drone
COLORE_FILO = [0.35, 0.35, 0.35]    # filo a piombo dal drone al pavimento: senza ombre è l'unico indizio della quota
COLORE_SFONDO = [0.82, 0.86, 0.91]  # dove non c'è nulla (cielo)

NOMI_CONTROLLORI = {'BS': 'backstepping', 'DFL': 'linearizzazione dinamica', 'GEO': 'controllo geometrico'}


def serve_ffmpeg():
    # Percorso di ffmpeg; se manca, errore subito e non a metà campagna
    percorso = shutil.which('ffmpeg')
    if percorso is None:
        raise RuntimeError("ffmpeg non trovato: serve per --video (winget install ffmpeg)")
    return percorso


def inquadratura(punti):
    # Matrice 3x4 della telecamera che vede tutti i punti da AZIMUT, ELEVAZIONE, puntata al centro della zona
    punti = np.array(punti, dtype=float)
    c = (punti.min(axis=0) + punti.max(axis=0))/2
    u = np.array([np.cos(ELEVAZIONE)*np.cos(AZIMUT), np.cos(ELEVAZIONE)*np.sin(AZIMUT), np.sin(ELEVAZIONE)])  # dal centro alla telecamera

    # La telecamera guarda lungo il suo asse z, con l'asse y verso l'alto dell'immagine
    z_c = -u
    y_c = np.array([0.0, 0.0, 1.0]) - z_c[2]*z_c
    y_c = y_c/np.linalg.norm(y_c)
    x_c = np.cross(y_c, z_c)

    # Distanza minima perché ogni punto (più il margine) stia dentro l'apertura, in orizzontale e in verticale
    tan_x = np.tan(ANGOLO_VISTA/2)
    tan_y = tan_x*ALTEZZA/LARGHEZZA
    q = punti - c
    profondita = q @ z_c
    d = max(DISTANZA_MIN,
            np.max((np.abs(q @ x_c) + MARGINE)/tan_x - profondita),
            np.max((np.abs(q @ y_c) + MARGINE)/tan_y - profondita))

    return np.column_stack([x_c, y_c, z_c, c + d*u])


def didascalia(d:dict):
    # Riferimento e scenario della prova, in una riga
    rif = d['riferimento']
    if rif['tipo'] == 'hover':
        eta0 = np.degrees(d['posa']['eta'])
        testo = "hover in "+str(np.round(rif['p'], 1).tolist())+" m, parte da "+str(np.round(d['posa']['p'], 1).tolist())+" m"
        if np.any(np.abs(eta0) > 1e-6):
            testo += ", eta "+str(np.round(eta0).astype(int).tolist())+"°"
    elif rif['tipo'] == 'otto':
        testo = "otto, ω = "+str(rif['omega'])+" rad/s"
    else:
        testo = "waypoint, velocità media "+str(rif['v_media'])+" m/s"

    if d['scenario'] == 'A':
        testo += "  |  vento "+str(d['aria']['vento'][0])+" m/s, raffica +"+str(d['aria']['raffica'][0])+" m/s a "+f"{d['aria']['t0']:.1f}"+" s"
    elif d['scenario'] == 'B':
        testo += "  |  carico ×"+str(d['fattore_carico'])+" (massa e inerzia)"
    elif d['scenario'] == 'C':
        testo += "  |  motori lenti, τm = "+str(int(round(1000*d['tau_m'])))+" ms"
    else:
        testo += "  |  nominale"
    return testo


def carattere(dimensione:int, grassetto:bool = False):
    for nome in (['segoeuib.ttf', 'arialbd.ttf'] if grassetto else ['segoeui.ttf', 'arial.ttf']):
        try:
            return ImageFont.truetype(nome, dimensione)
        except OSError:
            pass
    return ImageFont.load_default(dimensione)


class Ripresa:

    def __init__(self, sim, file:Path, descr:dict, rif, controllore:str, aria = None):
        # Da chiamare a scena pronta (drone nella posa iniziale), prima di sim.startSimulation()
        self.sim = sim
        self.file = Path(file)
        self.aria = aria if (aria is not None and aria.attiva) else None
        self.passi = max(1, int(round(1/(FPS*scenari.DT))))     # passi di simulazione per fotogramma
        self.fps = 1/(self.passi*scenari.DT)                    # [fotogrammi/s] esatti, così il video è a velocità reale

        # Riferimento campionato su tutta la prova: si disegna, e con la posa iniziale decide l'inquadratura
        tempi = np.arange(0.0, descr['durata'] + 1e-9, 0.02)
        punti_rif = np.array([rif(t)['p'] for t in tempi])
        punti = np.vstack([punti_rif, descr['posa']['p']])

        # Telecamera (sensore di visione) gestita da noi: si aggiorna solo quando serve un fotogramma
        opzioni = 1 + 2 + 4 + 128   # gestione esplicita, prospettiva, senza volume di vista disegnato, colore di sfondo nostro
        self.camera = sim.createVisionSensor(opzioni, [LARGHEZZA, ALTEZZA, 0, 0],
                                             [0.05, 100.0, float(ANGOLO_VISTA), 0.05, 0, 0, *COLORE_SFONDO, 0, 0])
        sim.setObjectAlias(self.camera, 'camera_video')
        sim.setObjectMatrix(self.camera, inquadratura(punti).flatten().tolist(), sim.handle_world)

        # I disegni fatti dal client NON spariscono con sim.loadScene: si tolgono quelli rimasti da prove precedenti
        # (interrotte, o di prima di questa correzione), e ogni ripresa toglie i suoi alla fine (togli_disegni)
        try:
            sim.removeDrawingObject(sim.handle_all)
        except Exception:
            pass

        # Riferimento intero (segmenti tra campioni consecutivi), scia del drone, punto desiderato, filo a piombo.
        # drawing_painttag: senza, i disegni si vedono nella finestra di CoppeliaSim ma non nella telecamera
        self.disegni = []
        segmenti = np.hstack([punti_rif[:-1], punti_rif[1:]])
        if len(segmenti):
            d_rif = sim.addDrawingObject(sim.drawing_lines + sim.drawing_painttag, 2, 0, -1, len(segmenti), COLORE_RIF)
            sim.addDrawingObjectItem(d_rif + sim.handleflag_addmultiple, segmenti.flatten().tolist())
            self.disegni.append(d_rif)
        n_fotogrammi = int(descr['durata']*self.fps) + 10
        self.d_scia = sim.addDrawingObject(sim.drawing_lines + sim.drawing_painttag, 3, 0, -1, n_fotogrammi, COLORE_SCIA)
        self.d_p_d = sim.addDrawingObject(sim.drawing_spherepts + sim.drawing_painttag, 0.06, 0, -1, 1, COLORE_RIF)
        self.d_filo = sim.addDrawingObject(sim.drawing_lines + sim.drawing_painttag, 1, 0, -1, 1, COLORE_FILO)
        self.disegni += [self.d_scia, self.d_p_d, self.d_filo]

        # Scritte
        self.titolo = descr['id']+"  ·  "+NOMI_CONTROLLORI.get(controllore, controllore)
        self.sottotitolo = didascalia(descr)
        self.c_titolo = carattere(30, True)
        self.c_testo = carattere(22)

        # ffmpeg riceve i fotogrammi RGB grezzi dallo standard input
        self.file.parent.mkdir(parents=True, exist_ok=True)
        comando = [serve_ffmpeg(), '-y', '-loglevel', 'error',
                   '-f', 'rawvideo', '-pix_fmt', 'rgb24', '-s', str(LARGHEZZA)+'x'+str(ALTEZZA), '-r', str(self.fps), '-i', '-',
                   '-c:v', 'libx264', '-crf', '20', '-pix_fmt', 'yuv420p', '-movflags', '+faststart', str(self.file)]
        self.ffmpeg = subprocess.Popen(comando, stdin=subprocess.PIPE, stderr=subprocess.PIPE)

        self.p_prec = None      # posizione all'ultimo fotogramma, per la scia
        self.ultima = None      # ultima immagine della telecamera, senza scritte
        self.stato = None       # (t, |e|, vento) dell'ultimo fotogramma

        return

    def fotogramma(self, k:int, t:float, p, p_d, ultimo:bool = False):
        # Chiamata a ogni passo k; riprende solo ogni self.passi passi, e all'ultimo passo della prova
        if k % self.passi != 0 and not ultimo:
            return
        sim = self.sim
        p = np.array(p, dtype=float)
        p_d = np.array(p_d, dtype=float)

        if self.p_prec is not None:
            sim.addDrawingObjectItem(self.d_scia, [*self.p_prec, *p])
        self.p_prec = p
        sim.addDrawingObjectItem(self.d_p_d, None)
        sim.addDrawingObjectItem(self.d_p_d, p_d.tolist())
        sim.addDrawingObjectItem(self.d_filo, None)
        sim.addDrawingObjectItem(self.d_filo, [*p, p[0], p[1], 0.0])

        sim.handleVisionSensor(self.camera)
        dati, risoluzione = sim.getVisionSensorImg(self.camera)
        immagine = np.frombuffer(dati, dtype=np.uint8).reshape(risoluzione[1], risoluzione[0], 3)
        self.ultima = immagine[::-1]        # OpenGL dà le righe dal basso verso l'alto

        vento = None if self.aria is None else float(np.linalg.norm(self.aria.velocita_vento(t)))
        self.stato = (t, float(np.linalg.norm(p_d - p)), vento)
        self.scrivi(self.componi())

        return

    def componi(self, esito:str = None, colore_esito = (255, 255, 255)):
        # Ultima immagine con titolo, didascalia, tempo, errore, vento ed eventualmente l'esito
        base = Image.fromarray(self.ultima).convert('RGBA')
        strato = Image.new('RGBA', base.size, (0, 0, 0, 0))
        disegno = ImageDraw.Draw(strato)
        t, e, vento = self.stato

        # In alto: chi e che cosa
        disegno.rounded_rectangle([12, 12, LARGHEZZA - 12, 96], radius=10, fill=(0, 0, 0, 150))
        disegno.text((28, 18), self.titolo, font=self.c_titolo, fill=(255, 255, 255, 255))
        disegno.text((28, 60), self.sottotitolo, font=self.c_testo, fill=(220, 220, 220, 255))

        # In basso a sinistra: tempo, errore di posizione, vento
        righe = ["t = "+f"{t:5.2f}"+" s", "|e| = "+f"{1000*e:6.1f}"+" mm"]
        if vento is not None:
            righe.append("vento "+f"{vento:4.1f}"+" m/s")
        alto = ALTEZZA - 22 - 32*len(righe)
        disegno.rounded_rectangle([12, alto - 10, 270, ALTEZZA - 12], radius=10, fill=(0, 0, 0, 150))
        for i, riga in enumerate(righe):
            disegno.text((28, alto + 32*i), riga, font=self.c_testo, fill=(255, 255, 255, 255))

        # In basso a destra: legenda
        x0 = LARGHEZZA - 250
        disegno.rounded_rectangle([x0, ALTEZZA - 86, LARGHEZZA - 12, ALTEZZA - 12], radius=10, fill=(0, 0, 0, 150))
        for i, (testo, colore) in enumerate([("riferimento", COLORE_RIF), ("drone", COLORE_SCIA)]):
            y = ALTEZZA - 64 + 32*i
            disegno.line([x0 + 16, y + 2, x0 + 60, y + 2], fill=tuple(int(255*c) for c in colore) + (255,), width=4)
            disegno.text((x0 + 74, y - 13), testo, font=self.c_testo, fill=(255, 255, 255, 255))

        # Esito, sull'ultimo fotogramma tenuto fermo
        if esito is not None:
            larghezza = disegno.textlength(esito, font=self.c_titolo)
            x = (LARGHEZZA - larghezza)/2
            disegno.rounded_rectangle([x - 20, 112, x + larghezza + 20, 162], radius=10, fill=(0, 0, 0, 180))
            disegno.text((x, 116), esito, font=self.c_titolo, fill=colore_esito + (255,))

        return Image.alpha_composite(base, strato).convert('RGB')

    def scrivi(self, immagine):
        self.ffmpeg.stdin.write(immagine.tobytes())

        return

    def togli_disegni(self):
        # Altrimenti riferimento, scia e p_d di questa prova comparirebbero nel video della prova dopo
        for h in self.disegni:
            try:
                self.sim.removeDrawingObject(h)
            except Exception:
                pass
        self.disegni = []

        return

    def chiudi(self, esito:dict, riassunto:dict):
        # Ultimo fotogramma fermo per ATTESA_FINALE s con l'esito, poi ffmpeg chiude il file
        self.togli_disegni()
        if self.ultima is not None:
            if esito['riuscita']:
                testo, colore = "riuscita  ·  e medio "+f"{1000*riassunto['e_medio']:.0f}"+" mm, e max "+f"{1000*riassunto['e_max']:.0f}"+" mm", (140, 230, 140)
            else:
                testo, colore = "FALLITA a t = "+f"{esito['t_fallimento']:.2f}"+" s: "+esito['motivo'], (255, 120, 110)
            immagine = self.componi(testo, colore)
            for _ in range(int(round(ATTESA_FINALE*self.fps))):
                self.scrivi(immagine)

        self.ffmpeg.stdin.close()
        errori = self.ffmpeg.stderr.read().decode(errors='replace')
        if self.ffmpeg.wait() != 0:
            raise RuntimeError("ffmpeg non ha scritto "+str(self.file)+": "+errori)

        return

    def annulla(self):
        # Prova interrotta: niente video a metà
        self.togli_disegni()
        try:
            self.ffmpeg.stdin.close()
        except OSError:
            pass
        self.ffmpeg.kill()
        self.ffmpeg.wait()
        self.file.unlink(missing_ok=True)

        return
