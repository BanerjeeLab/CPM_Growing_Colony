import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])
import CompuCellSetup


sim, simthread = CompuCellSetup.getCoreSimulationObjects()
CompuCellSetup.initializeSimulationObjects( sim, simthread)
        
        
        
#Add Python steppables here
steppableRegistry = CompuCellSetup.getSteppableRegistry()
        
from HST_Steppables import Initializer
Collective_Division_Initializer_Instance = Initializer(sim,_frequency=1)
steppableRegistry.registerSteppable( Collective_Division_Initializer_Instance )

from HST_Steppables import Crowded_Growth_Steppable
GrowthSteppableLinearInstance = Crowded_Growth_Steppable(sim,_frequency=1)
steppableRegistry.registerSteppable( GrowthSteppableLinearInstance )

from HST_Steppables import Sizer_Timer_Mitosis_Steppable
MitosisSteppableAdderInstance = Sizer_Timer_Mitosis_Steppable(sim,_frequency=10)
steppableRegistry.registerSteppable( MitosisSteppableAdderInstance )

from HST_Steppables import Density_Death_Steppable
DeathSteppableInstance = Density_Death_Steppable(sim,_frequency=10)
steppableRegistry.registerSteppable( DeathSteppableInstance )

from HST_Steppables import Extrusion
extrusionInstance = Extrusion(sim,_frequency=10)
steppableRegistry.registerSteppable( extrusionInstance )
 
 
 
        
CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        