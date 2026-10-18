import time
from coppeliasim_zmqremoteapi_client import RemoteAPIClient



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
    
    while (t := sim.getSimulationTime()) < 10:
        pos = sim.getObjectPosition(drone,sim.handle_world)
        z = pos[2]
        print("Simulation time:"+str(round(t,2))+"[s] " + "Altezza Drone:"+str(round(z,2)))
        sim.addForceAndTorque(radice,[0,0,massa*g],[0,0,0])
        sim.step()  # triggers next simulation step
    
    sim.stopSimulation()