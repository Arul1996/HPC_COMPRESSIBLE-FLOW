import numpy as np 
import para

# Central-difference method

def dfx_c(f):
    # Derivative of array 'f' w.r.t. x

    return (np.roll(f,-1,axis=0) - np.roll(f,1,axis=0))/(2*para.dx)

def dfy_c(f):
    # Derivative of array 'f' w.r.t. y

    return (np.roll(f,-1,axis=(para.dim-2)) - np.roll(f,1,axis=(para.dim-2)))/(2*para.dy)

def dfz_c(f):
    # Derivative of array 'f' w.r.t. z

    return (np.roll(f,-1,axis=(para.dim-1)) - np.roll(f,1,axis=(para.dim-1)))/(2*para.dz)

def d2fx_c(f):
    # 2nd derivative of array 'f' w.r.t. x

    return (np.roll(f,-1,axis=0) - ( 2*f ) + np.roll(f,1,axis=0))/(para.dx**2)

def d2fy_c(f):
    # 2nd derivative of array 'f' w.r.t. y

    return (np.roll(f,-1,axis=(para.dim-2)) - ( 2*f ) + np.roll(f,1,axis=(para.dim-2)))/(para.dy**2)

def d2fz_c(f):
    # 2nd derivative of array 'f' w.r.t. z

    return (np.roll(f,-1,axis=(para.dim-1)) - ( 2*f ) + np.roll(f,1,axis=(para.dim-1)))/(para.dz**2)

def dfxy_c(f):
    # Derivative of array 'f' w.r.t. x and y (d_x * d_y)

    return ((np.roll(f,(-1,-1),axis=(0,1)) - np.roll(f,(-1,1),axis=(0,1)) - np.roll(f,(1,-1),axis=(0,1)) + np.roll(f,(1,1),axis=(0,1)))/(4*para.dx*para.dy))

def dfyz_c(f):
    # Derivative of array 'f' w.r.t. y and z (d_y * d_z)

    return ((np.roll(f,(-1,-1),axis=((para.dim-2),(para.dim-1))) - np.roll(f,(-1,1),axis=((para.dim-2),(para.dim-1))) - np.roll(f,(1,-1),axis=((para.dim-2),(para.dim-1))) + np.roll(f,(1,1),axis=((para.dim-2),(para.dim-1))))/(4*para.dz*para.dy))

def dfzx_c(f):
    # Derivative of array 'f' w.r.t. z and x (d_z * d_x)

    return ((np.roll(f,(-1,-1),axis=(2,0)) - np.roll(f,(-1,1),axis=(2,0)) - np.roll(f,(1,-1),axis=(2,0)) + np.roll(f,(1,1),axis=(2,0)))/(4*para.dx*para.dz))


# Forward-difference method

def dfx_f(f,dx):
    # Derivative of array 'f' w.r.t. x

    return (np.roll(f,-1,axis=0) - f)/(dx)

def dfy_f(f,dy):
    # Derivative of array 'f' w.r.t. y

    return (np.roll(f,-1,axis=(para.dim-2)) - f)/(dy)

def dfz_f(f,dz):
    # Derivative of array 'f' w.r.t. z

    return (np.roll(f,-1,axis=(para.dim-1)) - f)/(dz)
