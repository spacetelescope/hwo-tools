import numpy as np

name = 'earth'    # Planet name from database
star = 'sun'      # Planet host star
distance = 10.0   # Distance to system (pc)
Nez = 1.0         # Number of exzodis (zodis)
Rp = 1.0          # Radius of planet (Earth Radii)
a = 1.0           # Semi-major axis (AU)
alpha = 90.       # Phase angle (deg)
Phi = 1./np.pi    # Lambertian phase function
MzV = 23.0        # Zodiacal light surface brightness (mag/arcsec**2)
MezV = 22.0       # Exozodiacal light surface brightness (mag/arcsec**2)