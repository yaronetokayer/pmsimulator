ACCEL_METHODS = ['tsc', 'cic', 'ngp']

def calculate_accel_ngp(x, y, accel_field):
    '''
    Nearest Grid Point (NGP) scheme for calculating acceleration of each particle

    Inputs:
    x - the x grid index of the particle
    x - the y grid index of the particle
    accel_field - accel_field attribute of Grid object

    Returns:
    accel_x, accel_y
    '''
    
    return accel_field[0][y, x], accel_field[1][y, x]