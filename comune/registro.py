"""
Salvataggio e lettura del log di una prova, in <cartella out>/<id>/:

    metadati.json   com'era fatta la prova, scritto una volta (testo, si legge a occhio)
    serie.npz       una riga per passo di controllo, con le stesse chiavi in tutte le prove
    interne.npz     le grandezze proprie del controllore, diverse per BS, DFL e GEO

Il log si scrive a prova finita, tutto in una volta: prima in <id>.parziale, poi la cartella si rinomina.
Così una prova interrotta non sembra fatta.
"""
import json
import shutil
from pathlib import Path

import numpy as np


def _in_json(x):
    # numpy -> tipi di base, per json
    if isinstance(x, dict):
        return {str(k): _in_json(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_in_json(v) for v in x]
    if isinstance(x, np.ndarray):
        return _in_json(x.tolist())
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    return x


def fatta(cartella:Path):
    return (Path(cartella) / 'metadati.json').exists()


def salva(cartella:Path, metadati:dict, serie:dict, interne:dict):
    cartella = Path(cartella)
    parziale = cartella.with_name(cartella.name + '.parziale')
    if parziale.exists():
        shutil.rmtree(parziale)
    parziale.mkdir(parents=True)

    with open(parziale / 'metadati.json', 'w', encoding='utf-8') as f:
        json.dump(_in_json(metadati), f, indent=2, ensure_ascii=False)
    np.savez_compressed(parziale / 'serie.npz', **serie)
    np.savez_compressed(parziale / 'interne.npz', **interne)

    if cartella.exists():
        shutil.rmtree(cartella)     # si rifà una prova già fatta (--rifai)
    parziale.rename(cartella)

    return


def carica(cartella:Path):
    # (metadati, serie, interne) di una prova
    cartella = Path(cartella)
    with open(cartella / 'metadati.json', encoding='utf-8') as f:
        metadati = json.load(f)
    with np.load(cartella / 'serie.npz') as z:
        serie = {k: z[k] for k in z.files}
    with np.load(cartella / 'interne.npz') as z:
        interne = {k: z[k] for k in z.files}
    return metadati, serie, interne
