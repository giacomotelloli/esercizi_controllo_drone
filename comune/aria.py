import numpy as np


class Aria:
    """
    Scenario A: attrito sulla velocità relativa all'aria, vento costante e raffica 1-coseno (planning §5.3).

        F_aria = -c_a m g (v - v_vento(t))
        v_vento(t) = vento + raffica/2 (1 - cos(2 pi (t - t0)/T_g))   per t0 <= t <= t0 + T_g

    c_a è in unità di m g per 1 m/s di velocità relativa. Senza descrizione (scenari N, B, C) la forza è nulla.
    """

    def __init__(self, descr:dict, massa:float, g:float = 9.81):
        self.attiva = descr is not None
        if not self.attiva:
            return
        self.c = descr['c_a']*massa*g                   # [N/(m/s)]
        self.vento = np.array(descr['vento'], dtype=float)
        self.raffica = np.array(descr['raffica'], dtype=float)
        self.t0 = descr['t0']
        self.T_g = descr['T_g']

        return

    def velocita_vento(self, t:float):
        v = self.vento.copy()
        if self.t0 <= t <= self.t0 + self.T_g:
            v = v + self.raffica/2*(1 - np.cos(2*np.pi*(t - self.t0)/self.T_g))
        return v

    def forza(self, t:float, v):
        # Forza [N] in terna mondo sul drone che si muove con velocità v [m/s]
        if not self.attiva:
            return np.zeros(3)
        return -self.c*(np.array(v) - self.velocita_vento(t))
