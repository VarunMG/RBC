from RBC_helper import *
from steady_state_helper3 import *


def longRun(Ra,Pr,alpha,Nx,Nz,T):
    #conduct long run at low Ra  
    RaBase = RBC_Problem(2000,7,1.5585,70,100)
    RaBase.initialize()
    RaBase.solve_system(500,True,False,True)
    print(RaBase.calc_Nu())
    RaBase.saveToFile('Ra2000Pr7alpha1.5585Nx70Nz100_T500.npy')
    RaBase.plot()

def baseSS(Ra,Pr,alpha,Nx,Nz,T,guessFile):
    #use long run to find the steady state at low Ra
    RaBase_SS = RBC_Problem(110000,7,1.5585,70,100)
    RaBase_SS.initialize()
    uArr, vArr, bArr, phiArr, dt = open_fields('Ra110000.00000000001Pr7alpha1.5585Nx70Nz100_SS.npy')
    RaBase_SS_guess = arrsToStateVec(phiArr, bArr)
    result = findSteadyState(RaBase_SS, RaBase_SS_guess, 2.0, 1e-6, 20,True)
    RaBase_SS.saveToFile('Ra110000Pr7alpha1.5585Nx70Nz100_SS.npy')

#using steady state at low Ra, follow branch upwards
uArr, vArr, bArr, phiArr, dt = open_fields('Ra459497.0Pr7alpha1.5585Nx70Nz100_SS.npy')
#for steady state helper 3: follow_branch(Pr,alpha,Ra_start,num_steps,Ra_step, Nx, Nz, startingGuess, starting_dt, tol):
# inputs for follow branch are like: follow_branch(Pr,alpha,Ra_start,Ra_end,Ra_step, Nx, Nz, startingGuess, starting_dt, tol)
starting_SS_state = arrsToStateVec(phiArr, bArr)
startingGuess = starting_SS_state
starting_dt = dt
RaVals, NuVals, steady_states = follow_branch(7,1.5585,459497,50, 1.05, 70, 100, startingGuess, starting_dt, 1e-6)

#uArr, vArr, bArr, phiArr, dt = open_fields('Ra100000Pr100alpha2.9634303295710804_optimalState.npy')
#starting_SS_state = arrsToStateVec(phiArr, bArr)
#startingGuess = starting_SS_state
#starting_dt = dt
#alpha_Vals, Nu_Vals, alphaOpt, NuOpt = findOptimalAlpha(100000,100,256,128,3.0,0.1,startingGuess,dt,1e-6,True)
#print("--------")
#print(alpha_Vals,Nu_Vals,alphaOpt,NuOpt)


#RaVals = np.arange(27000,50500,500)
#optAlphas = []
#optNus = []
#for i in range(len(RaVals)):
#    Ra = RaVals[i]
#    fileName = '/grad/gudibanda/RBC/steady_states/Pr7/primary_box/Nx256Nz128/Ra'+str(Ra)+'Pr7alpha1.5585Nx256Nz128_SS.npy'
#    uArr, vArr, bArr, phiArr, dt = open_fields(fileName)
#    starting_SS_state = arrsToStateVec(phiArr, bArr)
#    startingGuess = starting_SS_state
#    starting_dt = dt
#    alpha_Vals, Nu_Vals, alphaOpt, NuOpt = findOptimalAlpha(Ra,7,256,128,1.5585,0.1,startingGuess,dt,1e-3,True)
#    print("--------")
#    print(alpha_Vals,Nu_Vals,alphaOpt,NuOpt)
#    optAlphas.append(alphaOpt)
#    optNus.append(NuOpt)
#print('FINAL RESULTS')
#print("Ra Values")
#print(RaVals)
#print("optimal alpha values")
#print(optAlphas)
#print("optimal Nu values")
#print(optNus)

#uArr, vArr, bArr, phiArr, dt = open_fields('RB1_steady_states/Pr7/primary_box/Ra2000Pr7alpha1.5585Nx128Nz64data.npy')
#Pr7steady1 = RBC_Problem(3000,7,1.5585,128,64,'RB1')
#Pr7steady1.initialize()
#Pr7guess1 = arrsToStateVec(phiArr,bArr)
#iters1 = steady_state_finder(Pr7steady1,Pr7guess1,2,1e-2,50,False)
#follow_branch()

#alpha_vals, Nu_Vals, steady_states = optimize_alpha()
#print(alpha_vals)
#print(Nu_Vals)

# IH_test= RBC_Problem(1000,100,1.5585,256,128,'IH1')
# IH_test.initialize()
# IH_test.solve_system(30)
# IH_test.saveToFile('IH_time_stepping/Ra1000Pr100alpha1.5585Nx256Nz128.npy')


# alpha = 1.99945
# RBCLong = RBC_Problem(1150,100,alpha,256,128,'IH1')
# RBCLong.initialize()
# RBCLong.solve_system(100,True,False,True)
#RBCLong.saveToFile('temp_outputs/IH1_Ra4000Pr7alpha0.3141Nx512Nz256_T50.npy')


#tVals = RBCLong.tVals
#NuVals = RBCLong.NuVals
#file_name = 'temp_outputs/NuVals_IH1_Ra4000Pr7alpha0.3141Nx512Nz256_T50.npy'
#with open(file_name,'wb') as outFile:
#    np.save(outFile,NuVals)
#    np.save(outFile,tVals)
