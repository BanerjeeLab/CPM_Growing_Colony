from PlayerPython import * 
import CompuCellSetup
from PySteppables import *
import CompuCell
from PySteppablesExamples import MitosisSteppableBase 

import numpy as np
import random as random
import csv as csv
import time as time
import math as maths
import sys


######################################################################
##             Fundamental Paramters of the Simulation              ##
######################################################################
global CELL_ATTR_TRACKER
CELL_ATTR_TRACKER = {}
CELL_ATTR_TRACKER['MITOSIS_COND_DICT'] =    {
                                            'sizer_thresholds':{}, 
                                            'timer_thresholds':{}, 
                                            'timer_start_mcs' :{} 
                                            } 
CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'] = {}
CELL_ATTR_TRACKER['COLOR_CELLS'] = {}        

X_MAX = 1500
Y_MAX = 1500
EDGE_BUFFER = 50

SAVE_DATA_TO_CSV_INTERVAL = 500
END_SIM_IF_NO_CELLS       = True

folder_name_modifier = -1

### CHANGE THE PATH TO YOUR DESIRED OUTPUT
DATA_OUT_PATH = 'C:\CompuCell3D\CC3D_Workspace\output\\test\\' + str( folder_name_modifier ) + '\\'



## If  No_Death  or  Apoptosis_Off  are not in the DATA_OUT_PATH, then cells
#      will not die due to local density.
CALCULATE_MECH_PROB_APO = True
if 'No_Death' in DATA_OUT_PATH:
    CALCULATE_MECH_PROB_APO = False

if 'Apoptosis_Off' in DATA_OUT_PATH:
    CALCULATE_MECH_PROB_APO = False
    

######################################################################
##                         SIM PARAMETERS                           ##
######################################################################
SIM_CONSTS  = {}

SIM_CONSTS['subconfluent average area']     = 3200 

SIM_CONSTS['simulation initial cell length'] = 25 

SIM_CONSTS['timer threshold']       = 225 
SIM_CONSTS['timer threshold noise'] =   0.025 * SIM_CONSTS['timer threshold']

SIM_CONSTS['sizer threshold']       = 1377
SIM_CONSTS['sizer threshold noise'] =    0.025 * SIM_CONSTS['sizer threshold']

SIM_CONSTS['contact inhibition'] = 0.1
SIM_CONSTS['elastic modulous']   = 1.
SIM_CONSTS['growth rate'] = 0.5 * SIM_CONSTS['subconfluent average area'] / SIM_CONSTS['timer threshold'] 

SIM_CONSTS['density p_apo,max']            =   0.0015
SIM_CONSTS['density death sensitivity']    = 145.2     * 0.33**(-2.)
SIM_CONSTS['density death susceptibility'] =   0.01001 * 0.33**( 2.)

SIM_CONSTS['extrusion threshold'] = 0.25

######################################################################
##########           Func Defs for saving data        ################
######################################################################
def make_new_csv_file(filename, col_names):
    try:
        with open(filename, 'w') as file_obj:
            writer_obj = csv.writer(file_obj)
            writer_obj.writerow( col_names )
            return False
    except IOError:
        print("I/O error in make new csv "+filename)
        return True

def update_csv_file(filename, data):
    try:
        with open(filename, 'a+') as file_obj:
            writer_obj = csv.writer(file_obj)
            for line in data:
                writer_obj.writerow( line )
    except IOError:
        print("Failed to update csv "+filename)

def kill_simulation(self, reason_for_murder):
    for i in range(10):
        print(reason_for_murder)
        time.sleep(1)
    self.stopSimulation()  
    






######################################################################
######################################################################
class Initializer(SteppableBasePy):
    def __init__(self, _simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)           
        
    def make_cell(self):
        x, y, s = int(X_MAX/2), int(Y_MAX/2), SIM_CONSTS['simulation initial cell length']    
        cell = self.potts.createCell()
        self.cellField[x:x+s-1, y:y+s-1, 0] = cell
        
        CELL_ATTR_TRACKER[str(cell.id)] = {}          
        CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE'] = True

        cell.targetVolume = cell.volume 
        cell.type = 1
        cell.lambdaVolume = SIM_CONSTS['elastic modulous']   
    
    def start(self): 
        self.make_cell()
        



######################################################################
######################################################################
class Crowded_Growth_Steppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        global Ave_Vol_Data
        global Ave_Vol_Data_Filename
        Ave_Vol_Data_Filename = 'Ave_Vol_Data.csv'        
        Ave_Vol_Col_Names = ['Time',  'Count',  'Ave Vol',  'Total Vol',  'Ave dA/dt',  'Total dA/dt']
        file_failed = make_new_csv_file(DATA_OUT_PATH + Ave_Vol_Data_Filename, Ave_Vol_Col_Names)   
        
        if file_failed:
            kill_simulation(self, 'Failed to make ' + Ave_Vol_Data_Filename + 'Killing Sim')
        Ave_Vol_Data = []
    
        
    def step(self,mcs):    
        global Ave_Vol_Data
    
        count           = 0 
        total_area      = 0 
        added_dATdt     = 0 
        ave_area        = 0 
        ave_added_dATdt = 0 
                            
        for cell in self.cellList:   
            if CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE']:
                
                
                G = random.normalvariate( SIM_CONSTS['growth rate'], 0.05*SIM_CONSTS['growth rate'] )
                if G < 0:
                    G = SIM_CONSTS['growth rate']
                k = SIM_CONSTS['contact inhibition']
                dATdt = G * np.exp( -k * (cell.volume - cell.targetVolume)**2 )
                cell.targetVolume += dATdt
                
                added_dATdt += dATdt
                count       += 1
                total_area  += cell.volume

        ## END SIMULATION IF ALL CELLS ARE DEAD
        if (count == 0) and END_SIM_IF_NO_CELLS:
            kill_simulation(self, 'There are no cells left in the simulation...  Killing Sim')
         
               
        Ave_Vol_Data.append( [ mcs,  count,  ave_area,  total_area,  ave_added_dATdt,  added_dATdt ] )        
        
        if mcs % SAVE_DATA_TO_CSV_INTERVAL:
            update_csv_file( DATA_OUT_PATH + Ave_Vol_Data_Filename, Ave_Vol_Data )
            Ave_Vol_Data = []
            
        


            
######################################################################
######################################################################
class Sizer_Timer_Mitosis_Steppable(MitosisSteppableBase):  
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator,_frequency)
        
    def get_distance_to_center(self, x, y):
        return np.sqrt( (X_MAX*.5 - x)**2 + (Y_MAX*.5 - y)**2 )
        
    def get_distance_between_points(self, x1, y1, x2, y2):
        return np.sqrt( (x1 - x2)**2 + (y1 - y2)**2 )
        
    def start(self):   
        ##Initialize data save files 
        global Cell_Heritage_Data 
        global Heritage_Data_Filename
        Heritage_Data_Filename = 'Cell_Heritage_Data__Timer_Sizer_Mitosis.csv'
        Cell_Heritage_Col_Names = ['Time', 
                                        'Mother ID', 'Mother Type', 'Mother Volume', 'Mother Target Volume', 
                                                     'Mother X'      , 'Mother Y'   , 'Mother R'     , 
                                                     'Mother Sizer Threshold'       , 'Mother Timer Threshold', 'Mother Start of Timer Phase' ,
                                                     
                                        'Parent ID', 'Parent Type', 'Parent Volume', 'Parent Target Volume', 
                                                     'Parent X',       'Parent Y'   , 'Parent R', 
                                                     'Parent Sizer Threshold'       , 'Parent Timer Threshold',  
                                                     
                                        'Child ID' , 'Child Type',  'Child Volume',  'Child Target Volume',  
                                                     'Child X',     'Child Y',       'Child R', 
                                                     'Child Sizer Threshold',        'Child Timer Threshold' 
                                    ]
        file_failed = make_new_csv_file(DATA_OUT_PATH + Heritage_Data_Filename, Cell_Heritage_Col_Names)      
        if file_failed:
            kill_simulation(self, 'Failed to make ' + Heritage_Data_Filename + 'Killing Sim')
        Cell_Heritage_Data = []            
        
        global Cell_Cycle_Data_Filename 
        global Cell_Cycle_Data
        global Cell_Cycle_Col_Names
        Cell_Cycle_Data_Filename = 'Cell_Cycle_Data.csv'
        Cell_Cycle_Col_Names = [ 'Cycle Time'  , 'Cycle Area'   , 'Cycle Movement',              'Cycle R', 
                                 'Birth mcs'   , 'Birth Area'   , 'Birth X'   , 'Birth Y'   ,    'Birth R',
                                 'Division mcs', 'Division Area', 'Division X', 'Division Y', 'Division R',
                                 'Sizer Threshold', 
                                 'Timer Threshold',  
                                 'Start of Timer Phase' ] 
                                 
        file_failed = make_new_csv_file(DATA_OUT_PATH + Cell_Cycle_Data_Filename, Cell_Cycle_Col_Names)      
        if file_failed:
            kill_simulation(self, 'Failed to make ' + Cell_Cycle_Data_Filename + 'Killing Sim') 
        Cell_Cycle_Data = []
        
        
        
        ## Set sizer-timer thresholds for each cell at start of sim
        ## Set initial cell's values in CELL_CYCLE_DATA_DICT for data saving
        for cell in self.cellList:
            
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(cell.id)] = { col:None for col in Cell_Cycle_Col_Names }
            
            
            ### This is way overly complicated... sorry future me... 
            has_sizer_threshold = str(cell.id) in CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds']
            has_timer_threshold = str(cell.id) in CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds']
            
            if has_sizer_threshold:
                sizer_threshold = CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str( cell.id )] 
            else:
                sizer_threshold = int( random.normalvariate(SIM_CONSTS['sizer threshold'], SIM_CONSTS['sizer threshold noise']) )
                CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str( cell.id )]   = sizer_threshold
                
            
            if has_timer_threshold:
                timer_threshold = CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str( cell.id )] 
            else:
                timer_threshold = int( random.normalvariate(SIM_CONSTS['timer threshold'], SIM_CONSTS['timer threshold noise']) )
                CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str( cell.id )]   = timer_threshold
            
            
            
        
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( cell.id )]['Sizer Threshold'] = sizer_threshold
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( cell.id )]['Timer Threshold'] = timer_threshold
            
            
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( cell.id ) ]['Birth mcs']  = 0  
            ################ Target Vol because the cell is a small square at mcs=0
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( cell.id ) ]['Birth Area'] = cell.targetVolume  
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( cell.id ) ]['Birth X']    = cell.xCOM
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( cell.id ) ]['Birth Y']    = cell.yCOM
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( cell.id ) ]['Birth R']    = self.get_distance_to_center( cell.xCOM, cell.yCOM )
                                 
        
        
    def step(self,mcs):      
        global Cell_Heritage_Data 
        global Cell_Cycle_Data         
        global Cell_Cycle_Col_Names
        
        global CELL_ATTR_TRACKER
        
        cells_to_divide = []  

        ## Check if Mitosis Conditions Met for each cell
        for cell in self.cellList:
            if CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE']:
                
                ## Logic for Sizer-Timer Mitosis
                if str( cell.id ) in  CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs']:
                    time_past = mcs - CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs'][str(cell.id)]
                    time_thresh = CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str(cell.id)]
                    if time_past >= time_thresh:
                        cells_to_divide.append( cell )
                        CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs'].pop( str(cell.id), None )
                else:
                    if cell.volume > CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str(cell.id)]:
                        CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs'][str(cell.id)]         = mcs
                        CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(cell.id)]['Start of Timer Phase'] = mcs
                    
        
        
        
        for cell in cells_to_divide:       
        
            mother_id         = cell.id            
            mother_type       = cell.type         
            mother_vol        = cell.volume
            mother_target_vol = cell.targetVolume
            mother_x          = cell.xCOM
            mother_y          = cell.yCOM
            mother_R          = self.get_distance_to_center( mother_x, mother_y )
            mother_sizer_threshold   = CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( mother_id )]['Sizer Threshold']
            mother_timer_threshold   = CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( mother_id )]['Timer Threshold']
            mother_start_timer_phase = CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( mother_id )]['Start of Timer Phase']
        
            self.divideCellAlongMinorAxis(cell) #No calling the "cell" object after this            
            
            ############################################################################
            ## NOTE: "parentCell" is a terrible name because it is a daughter cell.    
            ##       CC3D gives the mother cell id to the "parentCell" and generates 
            ##       a new unique cell id for the "childCell".
            ############################################################################
            self.parentCell.targetVolume /= 2.     
            self.cloneParent2Child()  
            
            
            if CELL_ATTR_TRACKER['COLOR_CELLS'] == {}:
                if self.parentCell.type in [4,8]:
                    self.parentCell.type -= 3
                else:
                    self.parentCell.type = 1 + int(self.parentCell.type)   
                    
                    
            self.childCell.type = self.parentCell.type
            parent_id = self.parentCell.id
            child_id  = self.childCell.id
            birth_target_vol = self.childCell.targetVolume    
            
            # only need to set child because   mother_id = parent_id
            CELL_ATTR_TRACKER[str(child_id)] = {}                       
            CELL_ATTR_TRACKER[str(child_id)]['IS_ALIVE'] = True
            
            if CELL_ATTR_TRACKER['COLOR_CELLS'] != {}:
                ## It's a stupid declearation, but it stops python from making this entry a pointer.
                CELL_ATTR_TRACKER['COLOR_CELLS'][str( child_id )] = [ x for x in CELL_ATTR_TRACKER['COLOR_CELLS'][str( parent_id )] ]
            
            ################################################################################
            # 'Cycle Time'  ,   'Cycle  Area',      'Cycle Movement'     ,    'Cycle R',   #
            #                                                                              #
            # 'Birth mcs'   ,    'Birth Area',    'Birth X',    'Birth Y',    'Birth R',   #
            #                                                                              #
            # 'Division mcs', 'Division Area', 'Division X', 'Division Y', 'Division R',   #
            #                                                                              #
            # 'Sizer Threshold',                                                           #
            # 'Timer Threshold',                                                           #
            # 'Start of Timer Phase'                                                       #
            ################################################################################
            if str(cell.id) in CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT']:            
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Division mcs']  = mcs
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Division Area'] = mother_vol
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Division X']    = mother_x
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Division Y']    = mother_y
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Division R']    = mother_R
                
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Cycle Time']     = mcs - CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Birth mcs']
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Cycle Area']     = mother_vol - CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Birth Area']
                
                birth_x, birth_y = CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Birth X'] , CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Birth Y']                
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Cycle Movement'] = self.get_distance_between_points( mother_x, mother_y, birth_x, birth_y ) 
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Cycle R']        = CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Division R'] - CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Birth R'] 
                
                ## Store cell's cycle data in a temp list
                c_cyc_data = [ CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)][col] for col in Cell_Cycle_Col_Names ]
                ## Check that there are no "None" values in temp list
                if None in c_cyc_data:
                    None_cols = [ col for col in Cell_Cycle_Col_Names if CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)][col] == None ]
                    err_msg = 'None in Cell Cycle Data entry columns: '
                    for k in None_cols:
                        err_msg += k + ', '
                    kill_simulation(self, err_msg) # End sim if there are "None" values
                else:
                    Cell_Cycle_Data.append( c_cyc_data ) # Save the cycle data if all is kosher
                
            ## Initialize Cycle data for daughter cells
            for divided_cell in [ self.parentCell, self.childCell ]:               
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(divided_cell.id)] = { col:None for col in Cell_Cycle_Col_Names }
            
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( divided_cell.id ) ]['Birth mcs']  = mcs
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( divided_cell.id ) ]['Birth Area'] = divided_cell.volume
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( divided_cell.id ) ]['Birth X']    = divided_cell.xCOM
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( divided_cell.id ) ]['Birth Y']    = divided_cell.yCOM
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( divided_cell.id ) ]['Birth R']    = self.get_distance_to_center( divided_cell.xCOM, divided_cell.yCOM )
           
                
                
                sizer_threshold = int( random.normalvariate(SIM_CONSTS['sizer threshold'], SIM_CONSTS['sizer threshold noise']) )
                timer_threshold = int( random.normalvariate(SIM_CONSTS['timer threshold'], SIM_CONSTS['timer threshold noise']) )
                    
                CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str( divided_cell.id )]   = sizer_threshold
                CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str( divided_cell.id )]   = timer_threshold
                
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( divided_cell.id )]['Sizer Threshold'] = sizer_threshold
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( divided_cell.id )]['Timer Threshold'] = timer_threshold        
            

            
            Cell_Heritage_Data.append( [mcs, 
                mother_id, mother_type,  mother_vol, mother_target_vol  ,          
                    mother_x   ,  mother_y  , mother_R, 
                    mother_sizer_threshold, mother_timer_threshold, mother_start_timer_phase,
                            
                parent_id, self.parentCell.type, self.parentCell.volume, self.parentCell.targetVolume, 
                    self.parentCell.xCOM, self.parentCell.yCOM, CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( parent_id )]['Birth R'],
                    CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str( parent_id )], CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str( parent_id )],
                            
                child_id, self.childCell.type,  self.childCell.volume,  self.childCell.targetVolume,  
                    self.childCell.xCOM,  self.childCell.yCOM, CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( child_id )]['Birth R'],
                    CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str( child_id )], CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str( child_id )]
                                        ] )
        
           
        
        ## Save data
        if mcs % SAVE_DATA_TO_CSV_INTERVAL:
            update_csv_file( DATA_OUT_PATH + Heritage_Data_Filename, Cell_Heritage_Data )            
            Cell_Heritage_Data = []
            
            update_csv_file( DATA_OUT_PATH + Cell_Cycle_Data_Filename, Cell_Cycle_Data )            
            Cell_Cycle_Data = []
            
        
 
class Density_Death_Steppable(SteppableBasePy):
    def __init__(self, _simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
        
    def start(self): 
        global DENSITY_DEATH_RUNNING
        DENSITY_DEATH_RUNNING = True        

        global Local_Density_Data
        global Local_Density_Data_Filename
        Local_Density_Data_Filename = 'Local_Density_Data.csv'        
        Loc_Den_Col_Names = ['Time', 'ID', 'Type'
                     'Local Density', 'Number of Neighbors'
                     'Volume', 'Target Volume', 'Stress', 'Pressure', 'Perimeter Length',
                     'X', 'Y', 'R' ]
        file_failed = make_new_csv_file(DATA_OUT_PATH + Local_Density_Data_Filename, Loc_Den_Col_Names)   
        if file_failed:
            kill_simulation(self, 'Failed to make ' + Local_Density_Data_Filename + 'Killing Sim')  
        Local_Density_Data = []        
        
        global P_Apo_Density_Data
        global P_Apo_Density_Data_Filename
        P_Apo_Density_Data = []
        P_Apo_Density_Data_Filename = 'P_Apo_Density_Data.csv'  

        P_Apo_Density_Col_Names = Loc_Den_Col_Names
        if CALCULATE_MECH_PROB_APO:      
            file_failed = make_new_csv_file(DATA_OUT_PATH + P_Apo_Density_Data_Filename, P_Apo_Density_Col_Names)   
            if file_failed:
                kill_simulation(self, 'Failed to make ' + P_Apo_Density_Data_Filename + 'Killing Sim')                
                 
        global Vol_Den_Count_Data
        global Vol_Den_Count_Data_Filename
        Vol_Den_Count_Data_Filename = 'Size_Volume_Count_Data.csv'        
        Col_Names = ['Time', 'Average Density', 'Average Volume', 'Count' ]
        file_failed = make_new_csv_file(DATA_OUT_PATH + Vol_Den_Count_Data_Filename, Col_Names)   
        if file_failed:
            kill_simulation(self, 'Failed to make ' + Vol_Den_Count_Data_Filename + 'Killing Sim')  
        Vol_Den_Count_Data = []
        
    def death_prob(self, local_den, papo_max, alpha, rho_half):
        return papo_max/(1 + np.exp(-alpha * (local_den - rho_half)))
        
    def get_distance_to_center(self, x, y):
        return np.sqrt( (X_MAX*.5 - x)**2 + (Y_MAX*.5 - y)**2 )

    def step(self, mcs):        
        global Local_Density_Data   
        global Local_Density_Data_Filename        
        global Vol_Den_Count_Data
        global Vol_Den_Count_Data_Filename
        
        global P_Apo_Density_Data
                    
        count          = 0 
        summed_density = 0 
        total_vol      = 0            
        
        for cell in self.cellList:            
            # This ignores cells that are already dying
            if CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE'] and cell.volume > 1: 
                
                loc_density = cell.volume**-1.  
                num_neighbors = 0
                for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):   
                    if neighbor:                        
                        num_neighbors += 1                        
                        if (neighbor.volume > 0) and CELL_ATTR_TRACKER[str(neighbor.id)]['IS_ALIVE']: 
                            loc_density += neighbor.volume**-1.
                
                stress = (cell.volume - cell.targetVolume)**2                
                pressure = - 2 * SIM_CONSTS['elastic modulous'] * (cell.volume - cell.targetVolume)   
                cell_R = self.get_distance_to_center( cell.xCOM, cell.yCOM )
                
                perimeter_length = []
                for p in self.getCellBoundaryPixelList( cell ) :
                    perimeter_length.append( p )                    
                perimeter_length = len( perimeter_length )
                
                
                local_den_data_list = [mcs, cell.id, cell.type,
                                        loc_density, num_neighbors,
                                        cell.volume, cell.targetVolume, stress, pressure, perimeter_length,
                                        cell.xCOM, cell.yCOM, cell_R ]
                Local_Density_Data.append( local_den_data_list )    
                       
                count          += 1.
                summed_density += loc_density
                total_vol      += cell.volume
                
                if CALCULATE_MECH_PROB_APO:
                    in_X_bool = (cell.xCOM > EDGE_BUFFER and cell.xCOM < X_MAX-EDGE_BUFFER)
                    in_Y_bool = (cell.yCOM > EDGE_BUFFER and cell.yCOM < Y_MAX-EDGE_BUFFER)  
                    if in_X_bool and in_Y_bool:     
                        pmax     = SIM_CONSTS['density p_apo,max']
                        alpha    = SIM_CONSTS['density death sensitivity']
                        rho_Half = SIM_CONSTS['density death susceptibility']
                        if self.death_prob(loc_density, pmax, alpha, rho_Half) > random.random():
                            CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE'] = False
                            cell.targetVolume = 0
                            cell.lambdaVolume = 5    
                                
                        if not( CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE'] ):                             
                                P_Apo_Density_Data.append( local_den_data_list )
                            
        
                            
        Vol_Den_Count_Data.append([ mcs,  summed_density/count,  total_vol/count,  count])
        
        if mcs % SAVE_DATA_TO_CSV_INTERVAL:            
            update_csv_file( DATA_OUT_PATH + Local_Density_Data_Filename, Local_Density_Data )     
            update_csv_file( DATA_OUT_PATH + Vol_Den_Count_Data_Filename, Vol_Den_Count_Data )    
            
            Local_Density_Data = []
            Vol_Den_Count_Data = []   

            if CALCULATE_MECH_PROB_APO:
                update_csv_file( DATA_OUT_PATH + P_Apo_Density_Data_Filename, P_Apo_Density_Data )
                P_Apo_Density_Data = []      


    

  


class Extrusion(SteppableBasePy):        
    def __init__(self,_simulator,_frequency=20):
        SteppableBasePy.__init__(self,_simulator,_frequency)       
        
    def start(self):           
        global Extrusion_Data
        global Extrusion_Data_Filename
        Extrusion_Data_Filename = 'Extrusion_Data.csv'        
        Extrusion_Col_Names = ['Time', 'ID', 'Cell Type', 'Volume', 'Average Volume', 'Stress']
        file_failed = make_new_csv_file(DATA_OUT_PATH + Extrusion_Data_Filename, Extrusion_Col_Names)   
        if file_failed:
            kill_simulation(self, 'Failed to make ' + Extrusion_Data_Filename + 'Killing Sim') 
        Extrusion_Data = []
        
        
    def get_average_area(self):
        count, tot_area = 0., 0.
        for cell in self.cellList:
            tot_area += cell.volume
            count    += 1             
            
        return tot_area / count
        
        
        
    def step(self,mcs):
        global Extrusion_Data
        
        ave_area = self.get_average_area()
        for cell in self.cellList:            
            if cell.volume < (ave_area * SIM_CONSTS['extrusion threshold']):
                stress = (cell.volume - cell.targetVolume)**2.
                Extrusion_Data.append( [mcs, cell.id, cell.type, cell.volume, ave_area, stress] )
                self.deleteCell( cell )   
                    
        if mcs % SAVE_DATA_TO_CSV_INTERVAL and Extrusion_Data != []:
            update_csv_file( DATA_OUT_PATH + Extrusion_Data_Filename, Extrusion_Data )               
            Extrusion_Data = []
                







class ColorCellsByGrowthRateSteppable(SteppableBasePy):
    def __init__(self, _simulator, _frequency=1):
        SteppableBasePy.__init__(self, _simulator, _frequency)            
            
    def start(self):
        CELL_ATTR_TRACKER['COLOR_CELLS'] = {'Not_empty':[]}
        
        
    def step(self, mcs):    
        
        for cell in self.cellList:
            if str(cell.id) not in CELL_ATTR_TRACKER['COLOR_CELLS']:
                CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)] = []
        
        for cell in self.cellList:
            CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)].append( cell.volume )
            
            #if CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE']:
            if len( CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)] ) >= 10:
                dAdt = [] 
                if mcs%10 == 0:
                    for i,area_at_next_time in enumerate( CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)][1:] ):
                        dAdt.append( area_at_next_time - CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)][i] )
                    
                    #GAAAHHHHHHHH F#%$*% integer division!!  added 1.0 to fix
                    ave_G_rate = 1.0 *sum(dAdt) / len(dAdt)
                
                    clr_scale = [ 0, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8 ]                    
                    no_change_in_type = True
                    
                    if ave_G_rate >= clr_scale[-1]:
                        cell.type = 11
                        no_change_in_type = False
                    elif ave_G_rate <= -clr_scale[-1]:
                        cell.type = 32
                        no_change_in_type = False
                    
                    if ave_G_rate <= 0:
                        abs_Gr = ave_G_rate * -1.
                    else:
                        abs_Gr = ave_G_rate * 1.
                    
                    for i,s in enumerate( clr_scale[1:] ):
                    
                        if (clr_scale[i] <= abs_Gr) and (abs_Gr < s) and no_change_in_type:
                            if ave_G_rate < 0:
                                cell.type = 22 + i
                                no_change_in_type = False
                            else:  # 
                                cell.type = 21 - i
                                no_change_in_type = False
                    
                    
                        
                        
                #outside of the if mcs%10 --- delete oldest area value
                CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)].pop(0)    
                        
                        
                    
        
      
class ColorCellsByCellCyclePhaseSteppable(SteppableBasePy):
    def __init__(self, _simulator, _frequency=1):
        SteppableBasePy.__init__(self, _simulator, _frequency)
               
    def start(self):
        CELL_ATTR_TRACKER['COLOR_CELLS'] = {'Not_empty':[]}
            
    def step(self, mcs):    
        
        if mcs%5 == 0:
            for cell in self.cellList:
                if cell.volume <= CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str(cell.id)]:
                    size_thresh = CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str(cell.id)]
                    
                    frac_G1_completed = float( abs(cell.volume - size_thresh*0.5) / (size_thresh*0.5) )
                    G1_phase_type = int( frac_G1_completed *4 ) 
                    
                    #cell.type = 20
                    if cell.type != (11 + G1_phase_type):
                        cell.type = 11 + G1_phase_type
                elif str(cell.id) in CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs']:
                    time_past = mcs - CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs'][str(cell.id)]
                    time_thresh = CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str(cell.id)]
                    
                    frac_timer_phase_completed = 4 * float( time_past ) / float( time_thresh ) 
                    timer_phase_type = int( frac_timer_phase_completed )
                    
                    new_type = 14 + timer_phase_type                
                    if cell.type != new_type:
                        cell.type = new_type
                else:
                    cell.type = 14
