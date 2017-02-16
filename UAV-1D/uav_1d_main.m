% This program computes the trajectory of a 1D UAV
% - uav_1d_eom.m contains the Equations of Motion
% - controller.m contains the control law
% - planned_trajectory.m is used to return the desired final state and the
% initial condition
% - sys_params.m contains the system constants


% Initialize state
x0 = planned_trajectory(0);

% Time domain [s]
t0 = 0;
tf = 20;
params = sys_params; % put the system parameters in this structure

% Solve the ODE
options = odeset('RelTol',1e-5,'AbsTol',[1e-5 1e-5]);
[t,x]=ode45(@(t,x) uav_1d_eom(t,x,params),[t0 tf],x0,options);

% Recalculate control
for i = 1:length(t)
    u(i) = controller(t(i), x(i,:)', [1;0], params);
end
    
%plot(t,u,'b');
plot(t,x(:,1),'g');