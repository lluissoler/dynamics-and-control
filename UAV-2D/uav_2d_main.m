% This program computes the trajectory of a 2D UAV
% - uav_2d_eom.m contains the Equations of Motion
% - controller.m contains the control law
% - parameters.m contains the system constants

% Planned trajectory (uncomment desired trajectory):
    % trajhandle = @traj_step;
    % trajhandle = @traj_sine;
    % trajhandle = @traj_line;
    trajhandle = @traj_diamond;

% Define controller
controlhandle = @controller;

% Initialize system's parameters
params = parameters;

% Initialize state
des_start = trajhandle(0);
x0 = [des_start.pos; 0; des_start.vel; 0];

% Time domain [s]
t0 = 0;
tf = 10;

% Solve the ODE
options = odeset('RelTol',1e-5,'AbsTol',[1e-5 1e-5 1e-5 1e-5 1e-5 1e-5]);
[t,x]=ode45(@(t,x) uav_2d_eom(t,x,controlhandle,trajhandle,params),[t0 tf],x0,options);

% Create animation
figure;
hold on;
axis(gca,'equal');
yRange = max(x(:,1)) - min(x(:,1));
zRange = max(x(:,2)) - min(x(:,2));
axis([min(x(:,1))-0.1*yRange,max(x(:,1))+0.1*yRange,min(x(:,2))-0.1*zRange,max(x(:,2))+0.1*zRange]);

for i=1:round(length(t)/500):length(t)
    
    plot(x(i,1),x(i,2),'.b','MarkerSize',2);
    drawBody = plot(x(i,1),x(i,2),'or','MarkerSize',5);
    drawArms = line([x(i,1)+params.arm_length*cos(x(i,3));x(i,1)-params.arm_length*cos(x(i,3))],[x(i,2)+params.arm_length*sin(x(i,3));x(i,2)-params.arm_length*sin(x(i,3))]);
    timeText = text(min(x(:,1))-0.1*yRange+0.1,max(x(:,2))+0.1,sprintf('Time: %0.1f',t(i)));
    
    pause(0.0001);
    if i < length(t) - round(length(t)/500)
        delete(timeText);
        delete(drawBody);
        delete(drawArms);
    end
        
end









