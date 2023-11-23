def update_animation(frame, simulator, density_method, accel_method, adaptive_time):
    '''
    Update function for the animation.

    Parameters:
    - frame (int): The current frame number (automatically provided by FuncAnimation).
    - simulator (PMSimulator): The simulator object.
    - density_method (str): Density assignment method.
    - accel_method (str): Acceleration calculation method.
    - adaptive_time (bool): If True, uses adaptive time steps.
    '''
    nonlocal time_elapsed

    # Calculate the time step
    if adaptive_time:
        dt = simulator.compute_adaptive_dt(first=frame == 0)
    else:
        dt = simulator.simulation_settings.dt 

    # Update the grid and particle positions
    simulator.grid.assign_density(simulator.species_list, method=density_method)
    simulator.grid.compute_potential()

    for particle_pop in simulator.species_list:
        particle_pop.calculate_accels(simulator.grid, method=accel_method)
        particle_pop.update(dt)
    
    # Update the elapsed time
    time_elapsed += dt
    simulator.integration_time += dt

    # Update the plot for animation
    ax.clear()
    simulator.plot_snapshot(ax=ax, animate=True)