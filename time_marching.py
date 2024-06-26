import numpy as np
import dedalus.public as d3
import logging
import matplotlib.pyplot as plt
from matplotlib import cm
logger = logging.getLogger(__name__)


#def calcNu(b_var):
#    temp_field = b_var.allgather_data('g')
#    vert_means = np.mean(temp_field.T,axis=1)
#    return 1/(8*np.max(vert_means))

def calcNu(b_var):
    TAvg = d3.Average(b_var,coords['x']).evaluate()
    dzdTAvg = d3.Differentiate(TAvg,coords['z']).evaluate()
    dzdTAvgVals = dzdTAvg.allgather_data('g')
    #logger.info("Tz averages:")
    #logger.info(dzdTAvgVals)
    return -1*dzdTAvgVals[0][0]

def writeNu(fileName,tVals,NuVals):
    try:
        with open(fileName,'wb') as NuData:
            np.save(NuData,tVals)
            np.save(NuData,NuVals)
        return 1
    except:
        return 0
    
def calcCorrelation(field1,field2):
    field1Arr = field1.allgather_data('g').ravel()
    field2Arr = field2.allgather_data('g').ravel()
    costheta = np.dot(field1Arr,field2Arr)/(np.linalg.norm(field1Arr)*np.linalg.norm(field2Arr))
    return costheta

def writeFields(fileName,time,b_var,u_var,v_var):
    b_var.change_scales(1)
    u_var.change_scales(1)
    v_var.change_scales(1)
    with open(fileName,'wb') as fluidData:
        np.save(fluidData,time)
        np.save(fluidData,b_var.allgather_data('g').T)
        np.save(fluidData,u_var.allgather_data('g').T)
        np.save(fluidData,v_var.allgather_data('g').T)
    return 1

# Parameters
#Lx, Lz = 4, 1
#Lz = 2
alpha = 1.5585
Nx, Nz = 60, 100
Rayleigh = 100000
Prandtl = 7
dealias = 3/2
stop_sim_time = 500
timestepper = d3.RK222
max_timestep = 0.125
dtype = np.float64

# Bases
coords = d3.CartesianCoordinates('x', 'z')
dist = d3.Distributor(coords, dtype=dtype)
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(-1*np.pi/alpha, np.pi/alpha), dealias=dealias)
zbasis = d3.ChebyshevT(coords['z'], size=Nz, bounds=(-1, 1), dealias=dealias)

# Fields
phi = dist.Field(name='phi', bases=(xbasis,zbasis))
u = dist.Field(name='u', bases=(xbasis,zbasis))
v = dist.Field(name='v', bases=(xbasis,zbasis))
b = dist.Field(name='b', bases=(xbasis,zbasis))


tau_v1 = dist.Field(name='tau_v1', bases=xbasis)
tau_v2 = dist.Field(name='tau_v2', bases=xbasis)

tau_phi1 = dist.Field(name='tau_phi1', bases=xbasis)
tau_phi2 = dist.Field(name='tau_phi2', bases=xbasis)

tau_b1 = dist.Field(name='tau_b1', bases=xbasis)
tau_b2 = dist.Field(name='tau_b2', bases=xbasis)

tau_u1 = dist.Field(name='tau_u1', bases=xbasis)
tau_u2 = dist.Field(name='tau_u2', bases=xbasis)

# Substitutions
kappa = 4*(Rayleigh * Prandtl)**(-1/2)
nu = 4*(Rayleigh / Prandtl)**(-1/2)
x, z = dist.local_grids(xbasis, zbasis)
ex, ez = coords.unit_vector_fields(dist)
lift_basis = zbasis.clone_with(a=1/2, b=1/2) # First derivative basis
lift = lambda A, n: d3.Lift(A, lift_basis, n)

grad_v = d3.grad(v) + ez*lift(tau_v1,-1)
grad_phi = d3.grad(phi) + ez*lift(tau_phi1,-1)
grad_b = d3.grad(b) + ez*lift(tau_b1,-1) # First-order reduction
dz = lambda A: d3.Differentiate(A, coords['z'])
dx = lambda A: d3.Differentiate(A, coords['x'])

# Problem
# First-order form: "div(f)" becomes "trace(grad_f)"
# First-order form: "lap(f)" becomes "div(grad_f)"
problem = d3.IVP([phi, u, v, b, tau_v1, tau_v2, tau_phi1, tau_phi2, tau_b1, tau_b2], namespace=locals())
problem.add_equation("div(grad_v) + lift(tau_v2,-1) - phi= 0")
problem.add_equation("dt(phi) - nu*div(grad_phi)-  dx(dx(b)) + lift(tau_phi2,-1) = -dx(u*phi - v*lap(u))  ")
problem.add_equation("dt(b) - kappa*div(grad_b) + lift(tau_b2,-1) = -u*dx(b)-v*dz(b)")
problem.add_equation("dx(u) + dz(v)+ lift(tau_v1,-1) = 0", condition='nx!=0')
problem.add_equation("u = 0", condition='nx==0')
problem.add_equation("b(z=1) = -1")
problem.add_equation("v(z=1) = 0")
problem.add_equation("b(z=-1) = 1")
problem.add_equation("v(z=-1) = 0")
problem.add_equation("dz(v)(z=1) = 0")
problem.add_equation("dz(v)(z=-1) = 0")


# problem = d3.IVP([phi, v, b, tau_v1, tau_v2, tau_phi1, tau_phi2, tau_b1, tau_b2], namespace=locals())
# problem.add_equation("dt(phi) - nu*div(grad_phi) -  dx(dx(b)) + lift(tau_phi2,-1) = 0 ")
# problem.add_equation("dt(b) - kappa*div(grad_b) + lift(tau_b2,-1) = 0")
# problem.add_equation("div(grad_v) + lift(tau_v2,-1) - phi = 0")
# problem.add_equation("b(z=1) = -1")
# problem.add_equation("v(z=1) = 0")
# problem.add_equation("b(z=-1) = 1")
# problem.add_equation("v(z=-1) = 0")
# problem.add_equation("dz(v)(z=1) = 0")
# problem.add_equation("dz(v)(z=-1) = 0")


# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# Initial conditions
#b.fill_random('g', seed=42, distribution='normal', scale=1e-3) # Random noise
#b['g'] *= (1-z) * (1+z) # Damp noise at walls
b['g'] += - z # Add linear background
b['g'] += np.cos(alpha*x)*np.cos(np.pi*z/2)

# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=0.25, max_writes=50)
snapshots.add_task(b, name='buoyancy')
#snapshots.add_task(-d3.div(d3.skew(u)), name='vorticity')

# CFL
CFL = d3.CFL(solver, initial_dt=max_timestep, cadence=10, safety=0.5, threshold=0.05,
             max_change=1.5, min_change=0.5, max_dt=max_timestep)
CFL.add_velocity(u*ex+v*ez)


#solver.print_subproblem_ranks(dt=max_timestep)


# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=10)
#flow.add_property(np.sqrt(d3.dot(u,u))/nu, name='Re')
flow.add_property(1 + b*v/kappa,name='Nu')


##volume of box
volume = ((2*np.pi)/alpha)*2


# Main loop
startup_iter = 10

tVals = []
NuVals = []
#NuFileName = 'NuData.npy'

genFileName = 'Ra'+str(Rayleigh)+'Pr'+str(Prandtl)+'alpha'+str(alpha)+'Nx'+str(Nx)+'Nz'+str(Nz)+'_T' + str(stop_sim_time)
auxDataFile = genFileName + '_auxData/'
NuFileName = auxDataFile + genFileName + '_NuData.npy'
vertMeanFileName = auxDataFile + genFileName + '_vertMeans.npy'
fluidDataFileName = genFileName + '_runOutput/fluidData'

try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()
        solver.step(timestep)
        flow_Nu = flow.volume_integral('Nu')/volume
        #tVals.append(solver.sim_time)
        #NuVals.append(flow_Nu)
        #writeNu(NuFileName,tVals,NuVals)
        ####we try to make all the sin terms 0 in the horizontal direction to enforce symmetry
        b['c'][1::2,:] = 0
        v['c'][1::2,:] = 0
        phi['c'][1::2,:] = 0
        if (solver.iteration-1) % 10 == 0:
            #max_Re = flow.max('Re')
            #logger.info('Iteration=%i, Time=%e, dt=%e, max Re=%f' %(solver.iteration, solver.sim_time, timestep, max_Re))
            logger.info('Iteration=%i, Time=%e, dt=%e, Nu=%f' %(solver.iteration, solver.sim_time, timestep, flow_Nu))
        if (solver.iteration-1) %10 == 0:
            fileName = fluidDataFileName + str(round(100000*solver.sim_time)/100000) + '.npy'
            write = writeFields(fileName,solver.sim_time,b,u,v)
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()
    
def plot(field):
    field.change_scales(1)
    X,Z = np.meshgrid(x.ravel(),z.ravel())
    plt.contourf(X,Z,field['g'].T,cmap='seismic')
    plt.colorbar()
