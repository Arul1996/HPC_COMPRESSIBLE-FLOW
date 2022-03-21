import numpy as np
import para

class Compressible:

    def __init__(self):
        
        # Real-space velocity vector components
        self.ux = []
        self.uy = []
        self.uz = []
       
        # Non-linear calculation variables in momentum equation
        self.nlinx = [] 
        self.nliny = [] 
        self.nlinz = []
       
        # Real-space scalar fields
        self.rho = []
        self.p = []
        
        # Non-linear calculation variable in continuity equation
        self.nlinrho = []
        
        # rho * velocity mix vector components in real-space
        self.rho_ux = []
        self.rho_uy = []
        self.rho_uz = []
        
        # Total energy density variable
        self.E = []

        pass
    
    
    def set_arrays(self):
        # Used for setting arrays initially

        if(para.dim == 3):
            # Real-space velocity vector components
            self.ux = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])
            self.uy = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])
            self.uz = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])
     
            # Non-linear calculation variables in momentum equation
            self.nlinx = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])
            self.nliny = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])
            self.nlinz = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])
       
            # Real-space scalar fields
            self.rho = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])
            self.p = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])
                
            # Non-linear calculation variable in continuity equation
            self.nlinrho = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])
        
            # rho * velocity mix vector components in real-space
            self.rho_ux = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])
            self.rho_uy = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])
            self.rho_uz = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])

            # Total energy density variable
            self.E = np.ones([para.Nx+1,para.Ny+1,para.Nz+1])

            pass
    
        elif(para.dim == 2):
            # Real-space velocity vector components
            self.uy = np.ones([para.Ny+1,para.Nz+1])
            self.uz = np.ones([para.Ny+1,para.Nz+1])
        
            # Non-linear calculation variables in momentum equation
            self.nliny = np.ones([para.Ny+1,para.Nz+1])
            self.nlinz = np.ones([para.Ny+1,para.Nz+1])
       
            # Real-space scalar fields
            self.rho = np.ones([para.Ny+1,para.Nz+1])
            self.p = np.ones([para.Ny+1,para.Nz+1])
                
            # Non-linear calculation variable in continuity equation
            self.nlinrho = np.ones([para.Ny+1,para.Nz+1])
        
            # rho * velocity mix vector components in real-space
            self.rho_uy = np.ones([para.Ny+1,para.Nz+1])
            self.rho_uz = np.ones([para.Ny+1,para.Nz+1])

            # Total energy density variable
            self.E = np.ones([para.Ny+1,para.Nz+1])

            pass
        
        elif(para.dim == 1):
            # Real-space velocity vector components
            self.uz = np.ones([para.Nz+1])
     
            # Non-linear calculation variables in momentum equation
            self.nlinz = np.ones([para.Nz+1])
       
            # Real-space scalar fields
            self.rho = np.ones([para.Nz+1])
            self.p = np.ones([para.Nz+1])
                
            # Non-linear calculation variable in continuity equation
            self.nlinrho = np.ones([para.Nz+1])
        
            # rho * velocity mix vector components in real-space
            self.rho_uz = np.ones([para.Nz+1])

            # Total energy density variable
            self.E = np.ones([para.Nz+1])

            pass
              
        # size calculations

        # total_size = (float)(8 * 9 * (Nx * Ny * Nz)/(1024*1024))
        # total_size += (float)(8 * 12 * 2 *(Nkx * Nky * Nkz)/(1024*1024))
    
        # print("Total size = ", total_size)

        pass
    
    # @profile
    def init_hydro(self):
        # Initial profiles of variables
        
        if(para.dim == 3):
            
            # self.rho = np.sin(X_mesh)

            self.rho = 2* np.sin(para.X_mesh)
            self.ux = np.zeros([para.Nx+1,para.Ny+1,para.Nz+1])
            self.uy = np.zeros([para.Nx+1,para.Ny+1,para.Nz+1])
            self.uz = np.zeros([para.Nx+1,para.Ny+1,para.Nz+1])
        
            self.rho_ux = self.rho * self.ux
            self.rho_uy = self.rho * self.uy
            self.rho_uz = self.rho * self.uz
        
        if(para.dim == 2):
            
            # self.rho = 2* np.sin(para.Y_mesh)
            self.rho = np.ones([para.Ny+1,para.Nz+1])

            self.uy = np.zeros([para.Ny+1,para.Nz+1])
            self.uz = np.zeros([para.Ny+1,para.Nz+1])
            
            self.rho_uy = self.rho * self.uy
            self.rho_uz = self.rho * self.uz

        if(para.dim == 1):
            
            # self.rho = np.sin(para.Z)
            self.rho = np.ones(para.Nz+1)

            self.uz = np.zeros(para.Nz+1)
        
            self.rho_uz = self.rho * self.uz
            
        pass

    pass
