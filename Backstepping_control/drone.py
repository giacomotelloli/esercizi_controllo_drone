import numpy as np

class Drone: 

    def __init__(self, sim, path:str, param_motor:dict):

        base = sim.getObject(path+'/base')
        radice = sim.getObject(path)

        # Oggetto del drone in coppelia 
        self.drone = base 
        self.radice= radice 
        self.sim = sim

        # Calcolo Massa del Drone
        forme = sim.getObjectsInTree(radice,sim.sceneobject_shape)
        self.massa = 0 
        for handle in forme:
            if sim.getObjectInt32Param(handle,sim.shapeintparam_static) == 0 :
                # L'oggetto è dinamico e posso sommare la sua massa 
                self.massa += sim.getShapeMass(handle)
                print("Nome Pezzo:"+str(sim.getObjectAlias(handle)) + " Massa:"+str(sim.getShapeMass(handle)))
    
        print("Massa Totale:"+str(self.massa))

        # Inerzia del drone completo (serve al controllo d'assetto)
        self.calcola_inerzia(forme)
        print("Baricentro in terna body:"+str(self.baricentro))
        print("Inerzia J [kg*m^2]:\n"+str(self.J))

        # Mappatura dei motori e delle loro caratteristiche
        self.eliche = []
        self.pos_eliche = []
        for i in range(4):
            self.eliche.append(sim.getObject(path+'/propeller['+str(i)+']/respondable'))

            # Posizioni rispetto al drone
            self.pos_eliche.append(sim.getObjectPosition(self.eliche[i],self.radice)) 

        self.pos_eliche = np.array(self.pos_eliche) # 4x3, posizioni delle eliche nella terna del drone [m]

        # Parametri dei motori (scelti da noi, non si leggono dalla scena)
        self.k_f = param_motor['k_f']               # coefficiente di spinta [N/(rad/s)^2]
        self.k_m = param_motor['k_m']               # coefficiente di coppia di drag [N*m/(rad/s)^2]
        self.Omega_max = param_motor['Omega_max']   # velocità massima dei motori [rad/s]
        self.tau_m = param_motor['tau_m']           # costante di tempo dei motori [s]
        self.segni = np.array(param_motor['segni']) # verso di rotazione di ogni elica (+1/-1)
        self.T_max = 4 * self.k_f * self.Omega_max**2 # spinta totale con tutti i motori al massimo [N]

        # Discretizzazione esatta della dinamica del motore sul passo di simulazione
        self.dt = sim.getSimulationTimeStep()
        self.alpha = 1 - np.exp(-self.dt / self.tau_m)

        # Stato dei motori 
        self.Omega = np.zeros(4)

        x_eliche = self.pos_eliche[:,0]
        y_eliche = self.pos_eliche[:,1]

        # Righe della matrice di allocazione (mixer)
        riga_T = self.k_f * np.ones(4) 
        riga_tau_x = self.k_f * y_eliche
        riga_tau_y = -self.k_f*x_eliche
        riga_tau_z = self.k_m * self.segni

        # Matrice M 
        self.M = np.array([riga_T,riga_tau_x,riga_tau_y,riga_tau_z])
        self.M_inv = np.linalg.inv(self.M)

        # Inizializzazione Variabili 
        self.p = np.zeros(3) # posizione in sistema inerziale [m]
        self.v = np.zeros(3) # velocità in sistema inerziale [m/sec]
        self.eta = np.zeros(3) # phi, theta , psi [rad] , angoli di eulero della matrice R
        self.omega = np.zeros(3) # p, q, r velocità di rotazione attorno agli assi del drone (terna body) [rad/sec]
        self.R = np.eye(3) # Matrice di rotazione 

        return

    def calcola_inerzia(self, forme):
        # Inerzia dell'insieme dei pezzi dinamici attorno al baricentro comune, negli assi della terna body (self.drone)
        masse = []
        baricentri = []  # baricentro di ogni pezzo in terna body [m]
        inerzie = []     # inerzia di ogni pezzo attorno al proprio baricentro, assi body [kg*m^2]
        for handle in forme:
            if self.sim.getObjectInt32Param(handle,self.sim.shapeintparam_static) == 0 :
                # Inerzia attorno al baricentro del pezzo, con gli assi della terna del pezzo
                I_pezzo, com = self.sim.getShapeInertia(handle)
                I_pezzo = np.array(I_pezzo).reshape(3,3)
                c_pezzo = np.array([com[3], com[7], com[11]])

                # Posa del pezzo rispetto alla terna body
                mat = np.array(self.sim.getObjectMatrix(handle, self.drone)).reshape(3,4)
                R_pezzo = mat[:,:3]

                masse.append(self.sim.getShapeMass(handle))
                baricentri.append(R_pezzo @ c_pezzo + mat[:,3])
                inerzie.append(R_pezzo @ I_pezzo @ R_pezzo.T)

        masse = np.array(masse)
        baricentri = np.array(baricentri)
        self.baricentro = masse @ baricentri / np.sum(masse)  # baricentro complessivo in terna body [m]

        # Steiner: porto l'inerzia di ogni pezzo nel baricentro complessivo
        self.J = np.zeros((3,3))
        for m_i, c_i, I_i in zip(masse, baricentri, inerzie):
            d = c_i - self.baricentro
            self.J += I_i + m_i * (np.dot(d,d)*np.eye(3) - np.outer(d,d))

        return

    def aggiorna_stato(self):
        # Matrice del drone 
        m = self.sim.getObjectMatrix(self.drone, self.sim.handle_world)
        m = np.array(m)
        m = m.reshape(3,4)

        R = m[:,:3]
        self.R = R
        self.p = m[:,3]

        # Estraggo gli angoli di eulero dalla matrice 
        theta = -np.arcsin(R[2,0])
        phi = np.atan2(R[2,1],R[2,2])
        psi = np.atan2(R[1,0],R[0,0])
        angoli = [phi,theta,psi]
        self.eta = np.array(angoli)
        
        # Velocità del baricentro calcolate dal motore fisico, entrambe in terna mondo
        v_mondo , omega_mondo = self.sim.getVelocity(self.radice)
        self.v = np.array(v_mondo)                  # velocità lineare nel mondo [m/sec]
        self.omega = R.T @ np.array(omega_mondo)    # p, q, r in terna body [rad/sec]

        return 

    def applica_comando(self,T:float , tau:list):
        # Mixer: da spinta totale [N] e coppie in terna body [N*m] alle velocità dei motori
        u = np.array([T, tau[0], tau[1], tau[2]])
        Omega2 = self.M_inv @ u             # quadrati delle velocità [(rad/sec)^2]
        Omega2 = np.maximum(Omega2, 0)      # le eliche non possono spingere al contrario
        Omega_cmd = np.sqrt(Omega2)         # [rad/sec]

        self.comanda_motori(Omega_cmd)

        return


    def comanda_motori(self, Omega_cmd):
        # Saturazione: i motori girano solo tra 0 e Omega_max
        Omega_cmd = np.clip(np.array(Omega_cmd), 0, self.Omega_max)

        # Dinamica del primo ordine: la velocità reale insegue quella comandata
        self.Omega = self.Omega + self.alpha * (Omega_cmd - self.Omega)

        # Asse di spinta delle eliche (asse z del drone) espresso nel mondo
        asse = self.R[:,2]

        for i in range(4):
            f = self.k_f * self.Omega[i]**2                     # spinta dell'elica [N]
            Q = self.segni[i] * self.k_m * self.Omega[i]**2     # coppia di reazione [N*m]
            self.sim.addForceAndTorque(self.eliche[i], (f*asse).tolist(), (Q*asse).tolist())

        return