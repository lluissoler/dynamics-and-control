function [xdot] = uav_1d_eom(t,x,params)

% Differential equation for the 1D UAV
% x = [z, z_dot]
% x_dot = [z_dot, z_dot_dot]
% z_dot_dot = u/m - g

x_des = planned_trajectory(t,params);
u_contr = controller(t, x, x_des, params); % get the control value from the controller
u_actual = min(max(params.u_min, u_contr), params.u_max); % control saturation

xdot = [x(2);
        u_actual/params.mass - params.gravity];

end
