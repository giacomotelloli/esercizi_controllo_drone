"""
Campagna del backstepping in CoppeliaSim: le prove della griglia (doc/Campagna delle simulazioni.html),
3 riferimenti x 4 scenari x 3 livelli = 36 prove. Ogni prova scrive out/<id>_BS/:
metadati.json, serie.npz, interne.npz (comune/registro.py).

Tutto il resto è comune ai tre controllori: catalogo delle prove (comune/scenari.py), guadagni
(comune/taratura.py), ciclo di simulazione (comune/esecuzione.py), scena (comune/drone_control.ttt).

Uso, con CoppeliaSim aperto (la scena la carica lo script):
    python Backstepping_control/simula.py                             tutte le 36 prove
    python Backstepping_control/simula.py --test T1 --scenario N      solo alcune (anche --livello 1 2)
    python Backstepping_control/simula.py --elenco                    elenca le prove senza simulare
    python Backstepping_control/simula.py --rifai                     rifà anche quelle già in out/
    python Backstepping_control/simula.py --test T2 --livello 2 --video  anche il video, in out/<id>_BS/video.mp4
"""
import sys
from pathlib import Path

CARTELLA = Path(__file__).resolve().parent
sys.path.insert(0, str(CARTELLA.parent))

from comune.esecuzione import esegui_campagna
from controllore_backstepping import ControlloreBackstepping


if __name__ == "__main__":
    esegui_campagna('BS', ControlloreBackstepping, CARTELLA / 'out')
