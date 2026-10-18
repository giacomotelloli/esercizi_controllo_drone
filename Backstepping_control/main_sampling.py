import time
import numpy as np 
from coppeliasim_zmqremoteapi_client import RemoteAPIClient
np.set_printoptions(precision=3,suppress=True)

def leggi_stato(sim, drone , radice):
    # Matrice del drone 
    m = sim.getObjectMatrix(drone, sim.handle_world)
    m = np.array(m)
    m = m.reshape(3,4)

    R = m[:,:3]
    p = m[:,3]

    # Estraggo gli angoli di eulero dalla matrice 
    theta = -np.arcsin(R[2,0])
    phi = np.atan2(R[2,1],R[2,2])
    psi = np.atan2(R[1,0],R[0,0])

    # Estraggo le velocità di rotazione del corpo rispetto agli assi del mondo 
    # velocità del corpo in m/sec e angoli roll , pitch , yaw in rad/sec 
    v_body , omega_body = sim.getVelocity(radice)

    print(f"Matrice di rotazione:{R} \n Vettore di posizione: {p}")
    return 


if __name__== "__main__":
    print('Program started')

    # Connessione 
    client = RemoteAPIClient()
    sim = client.require('sim')

    drone = sim.getObject('/Quadcopter/base')
    print("Handle del drone:"+str(drone))
    # Run a simulation in stepping mode:
    sim.setStepping(True)
    sim.startSimulation()

    # Precalcolo fisica del drone 
    radice = sim.getObject('/Quadcopter')
    forme = sim.getObjectsInTree(radice,sim.sceneobject_shape)
    massa = 0 
    g = 9.81
    for handle in forme:
        if sim.getObjectInt32Param(handle,sim.shapeintparam_static) == 0 :
            # L'oggetto è dinamico e posso sommare la sua massa 
            massa += sim.getShapeMass(handle)
            print("Nome Pezzo:"+str(sim.getObjectAlias(handle)) + " Massa:"+str(sim.getShapeMass(handle)))

    print("Massa Totale:"+str(massa))

    # Ciclo della simulazione sincronizzata 
    while (t := sim.getSimulationTime()) < 10:
        pos = sim.getObjectPosition(drone,sim.handle_world)
        z = pos[2]
        print("Simulation time:"+str(round(t,2))+"[s] " + "Altezza Drone:"+str(round(z,2)))
        # Sto aggiungendo una forza virtuale 
        sim.addForceAndTorque(radice,[0,0,massa*g],[0,0,0])
        leggi_stato(drone=drone,radice=radice,sim=sim)
        sim.step()  # triggers next simulation step
    
    sim.stopSimulation()