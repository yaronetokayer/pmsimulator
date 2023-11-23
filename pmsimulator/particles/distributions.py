import numpy as np
# from astropy.cosmology import FlatLambdaCDM
# from astropy import units as u

# cosmo = FlatLambdaCDM(H0=69, Om0=0.3)

ALLOWED_DISTRIBUTIONS = ['uniform', 'normal', 'allones', 'nfw']

def initialize_uniform(n, bounds=(0, 1)):
    return np.random.uniform(low=bounds[0], high=bounds[1], size=n)

def initialize_normal(n, mean=0, std=1):
    return np.random.normal(mean, std, n)

def initialize_all_ones(n, scale=1.0):
    return scale * np.ones(n)

# def initialize_nfw(n, center=(0,0), r_200=1, c=5, z=0, cosmo=cosmo):
#     '''
#     Initialize the x and y coordinates of a particle distribution as an NFW
#     This would have to be a method of the Partcles class, since it intializes both the x and y coordinates simultaneously.

#     '''
#     # Define the range of values for x and y
#     x_values = np.linspace(-3, 3, 100)
#     y_values = np.linspace(-3, 3, 100)

#     # Create a meshgrid for x and y
#     x_mesh, y_mesh = np.meshgrid(x_values, y_values)

#     # Calculate the 2D PDF values
#     pdf_values = sigma_nfw(center, (x_mesh.flatten(), y_mesh.flatten()), r_200=r_200, c=c, cosmo=cosmo)
#     pdf_values = pdf_values.reshape(x_mesh.shape)

#     # Calculate the cumulative distribution function (CDF)
#     cdf_values = np.cumsum(pdf_values.flatten()) / np.sum(pdf_values)

#     # Generate random samples from a uniform distribution
#     uniform_samples = np.random.rand(n)

#     # Use the inverse of the CDF to map uniform samples to samples from the desired distribution
#     inverse_cdf_samples = np.interp(uniform_samples, cdf_values, np.arange(cdf_values.size))

#     # Map the 1D samples to 2D samples using the meshgrid
#     sample_indices = np.unravel_index(inverse_cdf_samples.astype(int), pdf_values.shape)
#     x_samples = x_mesh[sample_indices]
#     y_samples = y_mesh[sample_indices]

#     return x_samples, y_samples

# def sigma_nfw(center, coords, r_200, c=5, z=0, cosmo=cosmo):
# """
# NFW 2D projected surface mass density at r, given the concentration parameter and r_200
# e.g., Wright and Brainerd (1999) Eq. 11

# Inputs:
# center - tuple (x_center, y_center) of the center of the NFW profile
# coords - tuple (x, y) of coordinates at which to compute the density
# c - concentration parameter
# r_200 - array-like, radius of the halo inside which the mass density is 200*rho_c
# z - redshift of the halo (default is 0)
# cosmo - astropy cosmology class instantiation.
#         Default is FlatLambdaCDM(H0=69, Om0=0.3)

# Returns:
# sigma_nfw - 2D projected surface mass density at coords
# """

#     d_c = ( 200 / 3 ) * ( c**3 / ( np.log(1 + c) - ( c / (1 + c) ) ) )
#     h = (cosmo.H(z) * u.GeV * 6.5821e-25 * u.s).value # Hubble parameter in natural units
#     rho_c = ( 3 * h**2 ) / ( 8 * np.pi )
#     r_s = r_200 / c

#     factor = 2 * r_s * d_c * rho_c

#     x_rel = coords[0] - center[0]
#     y_rel = coords[1] - center[1]
    
#     x = ( np.sqrt(x_rel**2 + y_rel**2) / r_s ) # x is dimensionless distance
    
#     # x = 1
#     # returns zero for all other values
#     sigma_1 = np.where(x == 1, 1 / 3, np.zeros(x.shape))
    
#     # x < 1
#     # returns zero for all other values
#     sigma_l = ( 
#         np.power(x**2 - 1, -1, out=np.zeros(x.shape), where=x < 1 ) 
#         * ( 1 - 2 * np.arctanh( np.sqrt(( 1 - x ) / ( 1 + x), out=np.zeros(x.shape), where=x < 1  ) ) 
#            / np.sqrt( 1 - x**2, out=np.ones(x.shape), where=x < 1 ) 
#           )
#     )
    
#     # x > 1
#     # returns zero for all other values
#     sigma_g = ( 
#         np.power(x**2 - 1, -1, out=np.zeros(x.shape), where=x > 1 ) 
#         * ( 1 - 2 * np.arctan( np.sqrt(( x - 1 ) / ( 1 + x), out=np.zeros(x.shape), where=x > 1  ) ) 
#            / np.sqrt( x**2 - 1, out=np.ones(x.shape), where=x > 1 ) 
#           )
#     )
    
#     # Convert from avg surface density to total enclosed mass
#     sigma_nfw = factor * ( sigma_l + sigma_1 + sigma_g )
    
#     return sigma_nfw