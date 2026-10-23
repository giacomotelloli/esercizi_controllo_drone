import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from coppeliasim_zmqremoteapi_client import RemoteAPIClient

# La classe Drone è comune a tutti gli esercizi: si trova in comune/drone.py
sys.path.append(str(Path(__file__).resolve().parents[1]))
from comune.drone import Drone
