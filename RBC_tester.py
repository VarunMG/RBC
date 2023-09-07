from RBC_helper import *


###############################################################################

                                        #########
                                        #TESTING#
                                        #########


########################
### testing restarts ###
########################

def test_restart():
    Nx = 128
    Nz = 64
    test_start = RBC_Problem(5000,100,1.5585,Nx,Nz,'RB1')
    test_start.initialize()
    test_start.solve_system(100)
    test_start.saveToFile('test_save')
    
    test1 = RBC_Problem(5000,100,1.5585,Nx,Nz,'RB1')
    test1.initialize()
    test1.loadFromFile(1,'test_save')
    test1.solve_system(100)
    
    test_full = RBC_Problem(5000,100,1.5585,Nx,Nz,'RB1')
    test_full.initialize()
    test_full.solve_system(200)
    
    test1.u.change_scales(1)
    test1.v.change_scales(1)
    test1.b.change_scales(1)
    test1.phi.change_scales(1)
    
    test_full.u.change_scales(1)
    test_full.v.change_scales(1)
    test_full.b.change_scales(1)
    test_full.phi.change_scales(1)
    
    
    uErr = np.max(abs(test1.u['g']-test_full.u['g']))
    vErr = np.max(abs(test1.v['g']-test_full.v['g']))
    bErr = np.max(abs(test1.b['g']-test_full.b['g']))
    phiErr = np.max(abs(test1.phi['g']-test_full.phi['g']))
    print('u error: ', uErr)
    print('v error: ', vErr)
    print('b error: ', bErr)
    print('phi error: ', phiErr)
    
def test_restart_phi():
    Nx = 128
    Nz = 64
    test_start = RBC_Problem(5000,100,1.5585,Nx,Nz,'RB1')
    test_start.initialize()
    test_start.solve_system(100)
    test_start.saveToFile('test_save')
    
    uArr, vArr, bArr, phiArr, dt = open_fields('test_savedata.npy')
    
    test1 = RBC_Problem(5000,100,1.5585,Nx,Nz,'RB1')
    test1.initialize()
    test1.phi.load_from_global_grid_data(phiArr)
    test1.solve_system(100)
    
    test_full = RBC_Problem(5000,100,1.5585,Nx,Nz,'RB1')
    test_full.initialize()
    test_full.solve_system(200)
    
    test1.u.change_scales(1)
    test1.v.change_scales(1)
    test1.b.change_scales(1)
    test1.phi.change_scales(1)
    
    test_full.u.change_scales(1)
    test_full.v.change_scales(1)
    test_full.b.change_scales(1)
    test_full.phi.change_scales(1)
    
    
    uErr = np.max(abs(test1.u['g']-test_full.u['g']))
    vErr = np.max(abs(test1.v['g']-test_full.v['g']))
    bErr = np.max(abs(test1.b['g']-test_full.b['g']))
    phiErr = np.max(abs(test1.phi['g']-test_full.phi['g']))
    print('u error: ', uErr)
    print('v error: ', vErr)
    print('b error: ', bErr)
    print('phi error: ', phiErr)
    


#######################################
### testing the array manipulations ###
#######################################


def test_steady_finder():
    uArr1,vArr1,bArr1, phiArr1,dt1 = open_fields("test_files_curlcurl/steady_state_tests/Ra5000Pr100alpha1.5585Nx128Nz64T1000.npy")
    Nx = 128
    Nz = 64
    
    test1 = RBC_Problem(5000,100,1.5585,Nx,Nz,'RB1')
    test1.initialize()
    guess1 = arrsToStateVec(phiArr1, bArr1)
    iters1 = steady_state_finder(test1,guess1,2,1e-3,2,True)
    
    test1.b.change_scales(1)
    test1.phi.change_scales(1)
    
    steady1_b = test1.b.allgather_data('g')
    steady1_phi = test1.phi.allgather_data('g')
    
    
    error1_b = np.max(abs(steady1_b-bArr1))
    error1_phi = np.max(abs(steady1_phi-phiArr1))
    
    logging.info("number of iters: %i", iters1)
    logging.info("Error in b: %f",error1_b)
    logging.info("Error in phi: %f",error1_phi)
    logging.info("test 4 over \n ----------------------")
    
    
    uArr2,vArr2,bArr2, phiArr2,dt2 = open_fields("test_files_curlcurl/steady_state_tests/Ra5000Pr100alpha1.5585Nx128Nz64T250.npy")
    
    test2 = RBC_Problem(5000,100,1.5585,Nx,Nz,'RB1')
    test2.initialize()
    guess2 = arrsToStateVec(phiArr2, bArr2)
    iters2 = steady_state_finder(test2,guess2,2,1e-3,2,True)
    
    test2.b.change_scales(1)
    test2.phi.change_scales(1)
    
    steady2_b = test1.b.allgather_data('g')
    steady2_phi = test1.phi.allgather_data('g')
    
    error2_b = np.max(abs(steady2_b-bArr1))
    error2_phi = np.max(abs(steady2_phi-phiArr1))
    
    logging.info("number of iters: %i", iters2)
    logging.info("Error in b: %f",error2_b)
    logging.info("Error in phi: %f",error2_phi)
    logging.info("test 4 over \n ----------------------")
    
    return test2, bArr1,phiArr1