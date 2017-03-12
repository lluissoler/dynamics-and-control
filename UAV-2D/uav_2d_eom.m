function [xdot] = uav_2d_eom(t, x, controlhandle, trajhandle, params)

% Differential equation for the 2D UAV
% State:
%   [y; z; phi; y_dot; z_dot; phi_dot]

% Current state

current_state.pos = x(1:2);
current_state.rot = x(3);
current_state.vel = x(4:5);
current_state.omega = x(6);

% Desired state
desired_state = trajhandle(t, current_state);

% Control
[F, M] = controlhandle(t, current_state, desired_state, params);
u1 = 0.5*(F - M/params.arm_length);
u2 = 0.5*(F + M/params.arm_length);
u1_clamped = min(max(params.minF/2, u1), params.maxF/2);
u2_clamped = min(max(params.minF/2, u2), params.maxF/2);
F_clamped = u1_clamped + u2_clamped;
M_clamped = (u2_clamped - u1_clamped)*params.arm_length;
% F_clamped = 0;
% M_clamped = 0;

% ODE
xdot = [x(4);
        x(5);
        x(6);
        -F_clamped*sin(x(3))/params.mass;
        F_clamped*cos(x(3))/params.mass - params.gravity;
        M_clamped/params.Ixx];

    
end


