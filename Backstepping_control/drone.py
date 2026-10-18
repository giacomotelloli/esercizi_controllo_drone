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
        g = 9.81
        for handle in forme:
            if sim.getObjectInt32Param(handle,sim.shapeintparam_static) == 0 :
                # L'oggetto è dinamico e posso sommare la sua massa 
                self.massa += sim.getShapeMass(handle)
                print("Nome Pezzo:"+str(sim.getObjectAlias(handle)) + " Massa:"+str(sim.getShapeMass(handle)))
    
        print("Massa Totale:"+str(self.massa))

        # Inizializzazione Variabili 
        self.p = np.zeros(3) # posizione in sistema inerziale [m]
        self.v = np.zeros(3) # velocità in sistema inerziale [m/sec]
        self.eta = np.zeros(3) # phi, theta , psi [rad] , angoli di eulero della matrice R
        self.omega = np.zeros(3) # p, q, r velocità di rotazione rispetto agli assi del sistema inerziale
        self.R = np.eye(3) # Matrice di rotazione 


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
        
        # Estraggo le velocità di rotazione del corpo rispetto agli assi del mondo 
        # velocità del corpo in m/sec e angoli roll , pitch , yaw in rad/sec 
        v_body , omega_body = self.sim.getVelocity(self.radice)
        self.v = v_body
        self.omega = omega_body

        return 

    def applica_comando(self,T , tau):

        F_totale = T*self.R[:,2]
        return 
