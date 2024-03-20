from RBC_helper import *
from steady_state_helper3 import *

def saveArrs(uArr, vArr, bArr, phiArr,dt,outputFile):
    with open(outputFile,'wb') as outFile:
        np.save(outFile,uArr)
        np.save(outFile,vArr)
        np.save(outFile,bArr)
        np.save(outFile,phiArr)
        np.save(outFile,dt)
    return 0


def longRun(Ra,Pr,alpha,Nx,Nz,T):
    #conduct long run at low Ra  
    RaBase = RBC_Problem(Ra,Pr,alpha,Nx,Nz)
    RaBase.initialize()
    RaBase.solve_system(T,True,False,True)
    print(RaBase.calc_Nu())
    saveFile = 'Ra'+str(Ra)+'Pr'+str(Pr)+'alpha'+str(alpha)+'Nx'+str(Nx)+'Nz'+str(Nz)+'_T'+str(T)+'.npy'
    RaBase.saveToFile(saveFile)

def getSteady(Ra,Pr,alpha,Nx,Nz,T,guessFile):
    #use long run to find the steady state at low Ra
    RaBase_SS = RBC_Problem(Ra,Pr,alpha,Nx,Nz)
    RaBase_SS.initialize()
    uArr, vArr, bArr, phiArr, dt = open_fields(guessFile)
    RaBase_SS_guess = arrsToStateVec(phiArr, bArr)
    result = findSteadyState(RaBase_SS, RaBase_SS_guess, T, 1e-6, 20,True)
    saveFile = 'Ra'+str(Ra)+'Pr'+str(Pr)+'alpha'+str(alpha)+'Nx'+str(Nx)+'Nz'+str(Nz)+'_SS.npy'
    RaBase_SS.saveToFile(saveFile)

def refine(uArr,vArr,bArr,phiArr,alpha,Nx,Nz,newscale):
    ##run this using only ONE processor or else will not work
    print(uArr.shape)
    print("--------")
    coords = d3.CartesianCoordinates('x', 'z')
    dist = d3.Distributor(coords, dtype=np.float64)
    xbasis = d3.RealFourier(coords['x'], Nx, bounds=(-np.pi/alpha, np.pi/alpha), dealias=3/2)
    zbasis = d3.Chebyshev(coords['z'], Nz, bounds=(-1, 1), dealias=3/2)
    u = dist.Field(name='u', bases=(xbasis, zbasis))
    v = dist.Field(name='v', bases=(xbasis, zbasis))
    b = dist.Field(name='b', bases=(xbasis, zbasis))
    phi = dist.Field(name='phi', bases=(xbasis, zbasis))
    u.load_from_global_grid_data(uArr)
    v.load_from_global_grid_data(vArr)
    b.load_from_global_grid_data(bArr)
    phi.load_from_global_grid_data(phiArr)
    u.change_scales(newscale)
    v.change_scales(newscale)
    b.change_scales(newscale)
    phi.change_scales(newscale)
    newuArr = u['g']
    newvArr = v['g']
    newbArr = b['g']
    newphiArr = phi['g']
    print(newuArr.shape)
    return newuArr, newvArr, newbArr, newphiArr






###################
### refine grid ###
###################
'''
fileName = 'Ra90519.0Pr7alpha1.5585Nx420Nz600_SS.npy'
refinedFileName = 'Ra90519.0Pr7alpha1.5585Nx504Nz720_SS_refined.npy'
alpha = 1.5585

scaleFactor = 1.2
NxOld = 420
NzOld = 600


uArr, vArr, bArr, phiArr, dt = open_fields(fileName)
uRef, vRef, bRef, phiRef = refine(uArr,vArr,bArr, phiArr,alpha,NxOld,NzOld,scaleFactor)
saveArrs(uRef,vRef,bRef,phiRef,dt,refinedFileName)
'''
################
### Long run ###
################
#Ra=2000
#Pr=7
#alpha=1.5585
#Nx=70
#Nz=100
#T=500
#longRun(Ra,Pr,alpha,Nx,Nz,T)

#########################################
### For Finding a single steady state ###
#########################################
'''
Ra=90519
Pr=7
alpha=1.5585
Nx=504
Nz=720
T=2.0
guessFile = 'Ra90519.0Pr7alpha1.5585Nx504Nz720_SS_refined.npy'
getSteady(Ra,Pr,alpha,Nx,Nz,T,guessFile)

#using steady state at low Ra, follow branch upwards
#for steady state helper 3: follow_branch(Pr,alpha,Ra_start,num_steps,Ra_step, Nx, Nz, startingGuess, starting_dt, tol):
#inputs for follow branch are like: follow_branch(Pr,alpha,Ra_start,Ra_end,Ra_step, Nx, Nz, startingGuess, starting_dt, tol)
'''
startFile = 'Ra90519Pr7alpha1.5585Nx504Nz720_SS.npy'
uArr, vArr, bArr, phiArr, dt = open_fields(startFile)
starting_SS_state = arrsToStateVec(phiArr, bArr)
startingGuess = starting_SS_state
starting_dt = dt

startingRa = 90519
Pr = 7
alpha = 1.5585
Nx = 504
Nz = 720
numSteps = 15
RaStep = 1.1

RaVals, NuVals = follow_branch(Pr,alpha,startingRa,numSteps, RaStep, Nx, Nz, startingGuess, starting_dt, 1e-6)

#uArr, vArr, bArr, phiArr, dt = open_fields('Ra50000Pr7alpha1.5585Nx70Nz100_SS.npy')
#starting_SS_state = arrsToStateVec(phiArr, bArr)
#startingGuess = starting_SS_state
#starting_dt = dt
#alpha_Vals, Nu_Vals, alphaOpt, NuOpt = findOptimalAlpha(50000,7,70,100,1.5585,0.1,startingGuess,dt,1e-6,True)
#print("--------")
#print(alpha_Vals,Nu_Vals,alphaOpt,NuOpt)

