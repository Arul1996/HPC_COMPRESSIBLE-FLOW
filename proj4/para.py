import numpy as np
# from memory_profiler import profile

dim = 3                                     # Dimension of box
    
if(dim == 3):
    Nz = 32                                     # Number of grid points in z-direction
    Ny = 32                                     # Number of grid points in y-direction
    Nx = 32                                     # Number of grid points in x-direction

    Lx = np.pi                                  # Length of box in x-direction
    Ly = np.pi                                  # Length of box in y-direction
    Lz = np.pi                                  # Length of box in z-direction

    dx = Lx/Nx                                  # Length between two consecutive grid points in x-direction
    dy = Ly/Ny                                  # Length between two consecutive grid points in y-direction         
    dz = Lz/Nz                                  # Length between two consecutive grid points in z-direction

    X = np.arange(0,(Nx+1))*dx                  # Array consisting of points in x-direction
    Y = np.arange(0,(Ny+1))*dy                  # Array consisting of points in y-direction
    Z = np.arange(0,(Nz+1))*dz                  # Array consisting of points in z-direction

    X_mesh, Y_mesh, Z_mesh = np.meshgrid(X, Y, Z,indexing = 'ij')       # Meshgrids
    
if(dim == 2):
    Nz = 32                                     # Number of grid points in z-direction
    Ny = 32                                     # Number of grid points in y-direction
    
    Ly = np.pi                                  # Length of box in y-direction
    Lz = np.pi                                  # Length of box in z-direction

    dy = Ly/Ny                                  # Length between two consecutive grid points in y-direction         
    dz = Lz/Nz                                  # Length between two consecutive grid points in z-direction
    
    Y = np.arange(0,(Ny+1))*dy                  # Array consisting of points in y-direction
    Z = np.arange(0,(Nz+1))*dz                  # Array consisting of points in z-direction

    Y_mesh, Z_mesh = np.meshgrid(Y, Z,indexing = 'ij')                  # Meshgrids

if( dim == 1):
    Nz = 32                                     # Number of grid points in z-direction                    
    
    Lz = np.pi                                  # Length of box in z-direction

    dz = Lz/Nz                                  # Length between two consecutive grid points in z-direction

    Z = np.arange(0,(Nz+1))*dz                  # Array consisting of points in z-direction

tinit = 0                                   # Initial time
tfinal = 1                                  # Final time
dt = 0.00001                                # Single time step
tstep = np.array([0,0.5,1.572,2,3.14,4,5])  # Time step array for plotting

gamma = 5/3                                 # gamma = Cp/Cv
C = 1                                       # Equation of state constant, p = C * rho**gamma
mu = 0                                      # Dynamic viscosity

TempR = None                                # Temprary variables for calculations

n=1                                         # Used for rounding time in time-stepping
m=10 
while dt*m!=1:
    n=n+1
    m=m*10