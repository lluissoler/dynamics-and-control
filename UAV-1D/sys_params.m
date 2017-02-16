function [ params ] = sys_params()

% System constants
params.gravity = 9.8; % [m/s/s]
params.mass = 0.18; % [kg]
params.arm_length = 0.085; % [m]

% Thrust limits [N]
params.u_min = 0;
params.u_max = 2.5;

% Final state
params.z_des = 1;

end
