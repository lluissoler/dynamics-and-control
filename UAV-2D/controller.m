function [u1, u2] = controller(~, state, des_state, params)
%  CONTROLLER  Controller for the 2D quadrotor
%
%   Current state
%   state.pos = [y; z]
%   state.vel = [y_dot; z_dot]
%   state.rot = [phi]
%   state.omega = [phi_dot]
%
%   Desired state
%   des_state.pos = [y; z]
%   des_state.vel = [y_dot; z_dot]
%   des_state.acc = [y_ddot; z_ddot]
%

% PD Controller gains
kpz = 40;
kvz = 50;
kpphi = 20;
kvphi = 1;
kpy = 40;
kvy = 40;

% Trajectory controller
u11 = kvz*(des_state.vel(2) - state.vel(2));
u12 = kpz*(des_state.pos(2) - state.pos(2));
u1 = params.mass*(params.gravity + des_state.acc(2) + u11 + u12);

phi_com = -1/params.gravity*(des_state.acc(1) + kvy*(des_state.vel(1)-state.vel(1)) + kpy*(des_state.pos(1)-state.pos(1)));
u21 = kvphi*(-state.omega);
u22 = kpphi*(phi_com - state.rot);
u2 = (u21 + u22);

end
