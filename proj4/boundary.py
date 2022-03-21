import para
from compressible import Compressible

def boundary_dritchlet(compress = Compressible()):
    # At boundaries the value of rho, ux, uy and uz are zero
    
    if(para.dim == 3):
        
        compress.rho[0,:,:] = compress.rho[-1,:,:] = 0
        compress.rho[:,0,:] = compress.rho[:,-1,:] = 0
        compress.rho[:,:,0] = compress.rho[:,:,-1] = 0
    
        compress.ux[0,:,:] = compress.ux[-1,:,:] = 0
        compress.ux[:,0,:] = compress.ux[:,-1,:] = 0
        compress.ux[:,:,0] = compress.ux[:,:,-1] = 0
    
        compress.uy[0,:,:] = compress.uy[-1,:,:] = 0
        compress.uy[:,0,:] = compress.uy[:,-1,:] = 0
        compress.uy[:,:,0] = compress.uy[:,:,-1] = 0
    
        compress.uz[0,:,:] = compress.uz[-1,:,:] = 0
        compress.uz[:,0,:] = compress.uz[:,-1,:] = 0
        compress.uz[:,:,0] = compress.uz[:,:,-1] = 0

        pass
    
    if(para.dim == 2):
        
        compress.rho[0,:] = compress.rho[-1,:] = 0
        compress.rho[:,0] = compress.rho[:,-1] = 0
    
        compress.uy[0,:] = compress.uy[-1,:] = 0
        compress.uy[:,0] = compress.uy[:,-1] = 0
    
        compress.uz[0,:] = compress.uz[-1,:] = 0
        compress.uz[:,0] = compress.uz[:,-1] = 0

        pass
    
    if(para.dim == 1):
        
        compress.rho[0] = compress.rho[-1] = 0

        compress.uz[0] = compress.uz[-1] = 0
        
        pass
    
    pass
