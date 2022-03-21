from compressible import Compressible
from fns import time_advance_euler

compress = Compressible()           # Making of Object
compress.set_arrays()               # Setting of variables
compress.init_hydro()               # Initiating inital profiles at time tinit
time_advance_euler(compress)        # Time advance steps



