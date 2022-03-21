from compressible import Compressible
import para
from scipy.integrate import simpson as sn
from derivative import *
from boundary import boundary_dritchlet

def compute_RHS_rho(compress = Compressible()):
    # Computing scalar non-linear term in continuity equation

    if(para.dim == 3):
        compress.nlinrho[1:para.Nx,1:para.Ny,1:para.Nz] = dfx_c(compress.rho_ux)[1:para.Nx,1:para.Ny,1:para.Nz]
        compress.nlinrho[1:para.Nx,1:para.Ny,1:para.Nz] += dfy_c(compress.rho_uy)[1:para.Nx,1:para.Ny,1:para.Nz]
        compress.nlinrho[1:para.Nx,1:para.Ny,1:para.Nz] += dfz_c(compress.rho_uz)[1:para.Nx,1:para.Ny,1:para.Nz]    
        pass
    
    if(para.dim == 2):
        compress.nlinrho[1:para.Ny,1:para.Nz] += dfy_c(compress.rho_uy)[1:para.Ny,1:para.Nz]
        compress.nlinrho[1:para.Ny,1:para.Nz] += dfz_c(compress.rho_uz)[1:para.Ny,1:para.Nz]    
        pass
    
    if(para.dim == 1):
        compress.nlinrho[1:para.Nz] += dfz_c(compress.rho_uz)[1:para.Nz]    
        pass
    
    pass

# @profile
def compute_nlin_u(compress = Compressible()):
    # Computing vector non-linear term in momentum equation

    if(para.dim == 3):
        
        # Diagonal terms calculation
        compress.nlinz[1:para.Nx,1:para.Ny,1:para.Nz] = dfz_c(compress.rho_uz * compress.uz)[1:para.Nx,1:para.Ny,1:para.Nz]
        compress.nliny[1:para.Nx,1:para.Ny,1:para.Nz] = dfy_c(compress.rho_uy * compress.uy)[1:para.Nx,1:para.Ny,1:para.Nz]
        compress.nlinx[1:para.Nx,1:para.Ny,1:para.Nz] = dfx_c(compress.rho_ux * compress.ux)[1:para.Nx,1:para.Ny,1:para.Nz]
        
        # Off-diagonal terms calculation
        TempR = compress.rho_ux
        compress.rho_ux[1:para.Nx,1:para.Ny,1:para.Nz] = (compress.uy * compress.rho_ux)[1:para.Nx,1:para.Ny,1:para.Nz]
        compress.rho_uy[1:para.Nx,1:para.Ny,1:para.Nz] = (compress.uz * compress.rho_uy)[1:para.Nx,1:para.Ny,1:para.Nz]
        compress.rho_uz[1:para.Nx,1:para.Ny,1:para.Nz] = (TempR * compress.uz)[1:para.Nx,1:para.Ny,1:para.Nz]
    
        compress.nlinx[1:para.Nx,1:para.Ny,1:para.Nz] += dfy_c(compress.rho_ux)[1:para.Nx,1:para.Ny,1:para.Nz]
        compress.nlinx[1:para.Nx,1:para.Ny,1:para.Nz] += dfz_c(compress.rho_uz)[1:para.Nx,1:para.Ny,1:para.Nz]

        compress.nliny[1:para.Nx,1:para.Ny,1:para.Nz] += dfx_c(compress.rho_ux)[1:para.Nx,1:para.Ny,1:para.Nz]
        compress.nliny[1:para.Nx,1:para.Ny,1:para.Nz] += dfz_c(compress.rho_uy)[1:para.Nx,1:para.Ny,1:para.Nz]
    
        compress.nlinz[1:para.Nx,1:para.Ny,1:para.Nz] += dfx_c(compress.rho_uz)[1:para.Nx,1:para.Ny,1:para.Nz]
        compress.nlinz[1:para.Nx,1:para.Ny,1:para.Nz] += dfy_c(compress.rho_uy)[1:para.Nx,1:para.Ny,1:para.Nz] 
        
        pass
    
    if(para.dim == 2):
        
        # Diagonal terms calculation
        compress.nlinz[1:para.Ny,1:para.Nz] = dfz_c(compress.rho_uz * compress.uz)[1:para.Ny,1:para.Nz]
        compress.nliny[1:para.Ny,1:para.Nz] = dfy_c(compress.rho_uy * compress.uy)[1:para.Ny,1:para.Nz]
        
        # Off-diagonal terms calculation
        compress.rho_uy[1:para.Ny,1:para.Nz] = (compress.uz * compress.rho_uy)[1:para.Ny,1:para.Nz]

        compress.nliny[1:para.Ny,1:para.Nz] += dfz_c(compress.rho_uy)[1:para.Ny,1:para.Nz]
        compress.nlinz[1:para.Ny,1:para.Nz] += dfy_c(compress.rho_uy)[1:para.Ny,1:para.Nz] 
        
        pass
    
    if(para.dim == 1):
        
        # Diagonal terms calculation
        compress.nlinz[1:para.Nz] = dfz_c(compress.rho_uz * compress.uz)[1:para.Nz]
        
        pass
    
    pass

# @profile
def compute_p(compress = Compressible()):
    # Computing pressure from equation of state

    if(para.dim == 3):
        compress.p[1:para.Nx,1:para.Ny,1:para.Nz] = para.C * ((compress.rho[1:para.Nx,1:para.Ny,1:para.Nz])**para.gamma)
        pass
    
    if(para.dim == 2):
        compress.p[1:para.Ny,1:para.Nz] = para.C * ((compress.rho[1:para.Ny,1:para.Nz])**para.gamma)
        pass
    
    if(para.dim == 1):
        compress.p[1:para.Nz] = para.C * ((compress.rho[1:para.Nz])**para.gamma)
        pass
    
    pass

# @profile
def compute_RHS_rho_u(compress = Compressible()):
    # Computing R.H.S of momentum equation

    if( para.dim == 3):  

        compress.nlinx = - compress.nlinx - (dfx_c(compress.p)) 
        compress.nlinx += para.mu*((d2fx_c(compress.ux) + d2fy_c(compress.ux) + d2fz_c(compress.ux)) + ( (1/3) * ( d2fx_c(compress.ux) + dfxy_c(compress.uy) + dfzx_c(compress.uz))))
    
        compress.nliny = - compress.nliny - (dfy_c(compress.p)) 
        compress.nliny += para.mu*((d2fx_c(compress.uy) + d2fy_c(compress.uy) + d2fz_c(compress.uy)) + ( (1/3) * ( dfxy_c(compress.ux) + d2fy_c(compress.uy) + dfyz_c(compress.uz))))
    
        compress.nlinz = - compress.nlinz - (dfz_c(compress.p)) 
        compress.nlinz += para.mu*((d2fx_c(compress.uz) + d2fy_c(compress.uz) + d2fz_c(compress.uz)) + ( (1/3) * ( dfzx_c(compress.ux) + dfyz_c(compress.uy) + d2fz_c(compress.uz))))
    
        pass
    
    if(para.dim == 2):  

        compress.nliny = - compress.nliny - (dfy_c(compress.p)) 
        compress.nliny += para.mu*((d2fy_c(compress.uy) + d2fz_c(compress.uy)) + ( (1/3) * ( d2fy_c(compress.uy) + dfyz_c(compress.uz))))
    
        compress.nlinz = - compress.nlinz - (dfz_c(compress.p)) 
        compress.nlinz += para.mu*(( d2fy_c(compress.uz) + d2fz_c(compress.uz) ) + ( (1/3) * ( dfyz_c(compress.uy) + d2fz_c(compress.uz))))
        
        pass
    
    if(para.dim == 1):
        
        compress.nlinz = - compress.nlinz - (dfz_c(compress.p)) 
        compress.nlinz += para.mu*( d2fz_c(compress.uz)  + ( (1/3) * d2fz_c(compress.uz)))
        
        pass
    
    pass

# @profile
def time_advance_single_step(compress = Compressible()):
    
    compute_RHS_rho(compress)
    compute_nlin_u(compress)
    compute_p(compress)
    compute_RHS_rho_u(compress)
        
    if( para.dim == 3):
        
        # Computing rhok in time = tinit + dt 
        compress.rho += -compress.nlinrho * para.dt
    
        # Computing rho_u in time = tinit + dt 
        compress.rho_ux += compress.nlinx * para.dt
        compress.rho_uy += compress.nliny * para.dt
        compress.rho_uz += compress.nlinz * para.dt
        pass
    
    if( para.dim == 2):
        
        # Computing rhok in time = tinit + dt 
        compress.rho += -compress.nlinrho * para.dt
    
        # Computing rho_u in time = tinit + dt 
        compress.rho_uy += compress.nliny * para.dt
        compress.rho_uz += compress.nlinz * para.dt
        pass
    
    if( para.dim == 1):
        
        # Computing rhok in time = tinit + dt 
        compress.rho += -compress.nlinrho * para.dt
    
        # Computing rho_u in time = tinit + dt 
        compress.rho_uz += compress.nlinz * para.dt
        pass
    
    # Boundary condition
    boundary_dritchlet(compress)   
    
    pass

# @profile
def time_advance_euler(compress = Compressible()):
    
    t = para.tinit
   
    while(t <= para.tfinal):

        # Sequence is same as calculation
        time_advance_single_step(compress)
        
        if(para.dim == 3):
            # Computing respective variables
            compress.ux[1:para.Nx,1:para.Ny,1:para.Nz] = compress.rho_ux[1:para.Nx,1:para.Ny,1:para.Nz] / compress.rho[1:para.Nx,1:para.Ny,1:para.Nz]
            compress.uy[1:para.Nx,1:para.Ny,1:para.Nz] = compress.rho_uy[1:para.Nx,1:para.Ny,1:para.Nz] / compress.rho[1:para.Nx,1:para.Ny,1:para.Nz]
            compress.uz[1:para.Nx,1:para.Ny,1:para.Nz] = compress.rho_uz[1:para.Nx,1:para.Ny,1:para.Nz] / compress.rho[1:para.Nx,1:para.Ny,1:para.Nz]
            pass
        
        if(para.dim == 2):
            # Computing respective variables
            compress.uy[1:para.Ny,1:para.Nz] = compress.rho_uy[1:para.Ny,1:para.Nz] / compress.rho[1:para.Ny,1:para.Nz]
            compress.uz[1:para.Ny,1:para.Nz] = compress.rho_uz[1:para.Ny,1:para.Nz] / compress.rho[1:para.Ny,1:para.Nz]
            pass
        
        if(para.dim == 1):
            # Computing respective variables
            compress.uz[1:para.Nz] = compress.rho_uz[1:para.Nz] / compress.rho[1:para.Nz]
            pass
        
        # plt.plot(t,total_energy(compress))
        total_energy_density(compress)
        print(total_energy(compress))
        t=np.round(t+para.dt,para.n)         
        pass
    # plt.show()

    
def total_energy_density(compress = Compressible()):
    # Total energy density 
    
    if( para.dim == 3):
        compress.E = (para.C * ((compress.rho**(para.gamma-1))/(para.gamma-1))) + 1/2 * (compress.ux**2 + compress.uy**2 + compress.uz**2)
        
    if(para.dim == 2):
        compress.E = (para.C * ((compress.rho**(para.gamma-1))/(para.gamma-1))) + 1/2 * (compress.uy**2 + compress.uz**2)
    
    if(para.dim == 1):
        compress.E = (para.C * ((compress.rho**(para.gamma-1))/(para.gamma-1))) + (1/2 * (compress.uz**2))
    pass

def total_energy(compress = Compressible()):
    # Total energy integral
    
    if(para.dim == 1):
        return sn(compress.E,para.Z)
    
    elif(para.dim == 2):
        return sn(sn(compress.E,para.Z),para.Y)
    
    elif(para.dim == 3):
        return sn(sn(sn(compress.E,para.Z),para.Y),para.X)
