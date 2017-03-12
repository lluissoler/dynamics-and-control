function [params] = parameters()

% Returns a structure with the constants of the system

params.gravity = 9.81;
params.mass = 0.2;
params.Ixx = 0.00025;
params.arm_length = 0.1;

params.minF = 0;
params.maxF = 2*params.mass*params.gravity;

end
