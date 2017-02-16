This program computes the trajectory of a 1D UAV:
- Execute uav_1d_main to generate the trajectory. x(:,1) is the height (z) and x(:,2) the vertical velocity
- uav_1d_eom.m contains the Equations of Motion
- controller.m contains the control law
- planned_trajectory.m is used to return the desired final state and the initial condition
- sys_params.m contains the system constants
