
# def init_animation(ax, simulator):
#     '''
    
#     '''
#     # Draw initial condition
#     for particle_pop in self.species_list:
#             scatter_kwargs = {
#                 'marker': kwargs.get('marker', 'o'),
#                 's': kwargs.get('s', 2),
#                 'alpha': kwargs.get('alpha', 0.5),
#                 'label': particle_pop.species_name
#             }
#             ax.scatter(
#                 particle_pop.positions.x.data, particle_pop.positions.y.data, **scatter_kwargs
#             )

#     ax.set_title(f't={self.integration_time:.4f}')
#     ax.set_xlabel('X Position')
#     ax.set_ylabel('Y Position')
#     ax.set_xlim(0, self.simulation_settings.domain_size)
#     ax.set_ylim(0, self.simulation_settings.domain_size)

# def update_animation(frame, simulator, density_method, accel_method, adaptive_time):
#     '''
#     Update function for the animation.

#     Parameters:
#     - frame (int): The current frame number (automatically provided by FuncAnimation).
#     - simulator (PMSimulator): The simulator object.
#     - density_method (str): Density assignment method.
#     - accel_method (str): Acceleration calculation method.
#     - adaptive_time (bool): If True, uses adaptive time steps.
#     '''
#     nonlocal time_elapsed

#     # Calculate the time step
#     if adaptive_time:
#         dt = simulator.compute_adaptive_dt(first=frame == 0)
#     else:
#         dt = simulator.simulation_settings.dt 

#     # Update the grid and particle positions
#     simulator.grid.assign_density(simulator.species_list, method=density_method)
#     simulator.grid.compute_potential()

#     for particle_pop in simulator.species_list:
#         particle_pop.calculate_accels(simulator.grid, method=accel_method)
#         particle_pop.update(dt)
    
#     # Update the elapsed time
#     time_elapsed += dt
#     simulator.integration_time += dt

#     if frame is not None:
#         # Update the plot for animation
#         ax.clear()
#         simulator.plot_snapshot(ax=ax, animate=True)