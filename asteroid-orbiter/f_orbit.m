function xdot=f_orbit(t,x,R0,Rs,G,Msc,Msun,muast,i,W,P0,A,w,plbd,pphi,C20,C22)

rsun = [-Rs*cos(W*t); -Rs*sin(W*t)*cos(i); Rs*sin(W*t)*sin(i)]; % Position of the Sun
% Rzplbd= [cos(plbd) -sin(plbd) 0; sin(plbd) cos(plbd) 0; 0 0 1];
% Rypphi= [cos(pphi) 0 sin(pphi);0 1 0; -sin(pphi) 0 cos(pphi)];
% rsun = (Rzplbd*Rypphi)*rsun;

r = [x(1); x(2); x(3)]; % Position of the spacecraft

phi = pi/2-asin(x(3)/sqrt(norm(r)));
lbd = atan2(x(2),x(1));

dUdr = - 3/norm(r)^4*((3/2*sin(phi)^2-1/2)*C20 + 3*cos(phi)^2*C22*cos(2*lbd));
%dUdr = -1/norm(r)^2;
%dUdr = -1/norm(r)^2 - 3/norm(r)^4*((3/2*sin(phi)^2-1/2)*C20 + 3*cos(phi)^2*C22*cos(2*lbd));
dUdp = 1/norm(r)^4*(3*C20*cos(phi)*sin(phi)-6*C22*cos(2*lbd)*cos(phi)*sin(phi));
%dUdp = 0;
dUdl = -1/norm(r)^4*(6*C22*cos(phi)*sin(2*lbd));
%dUdl = 0;
rdd = muast*[dUdr; dUdp; dUdl]; % Acceleration vector in body fixed frame

%Rzlbd = [cos(lbd-w*t) -sin(lbd-w*t) 0; sin(lbd-w*t) cos(lbd-w*t) 0; 0 0 1];
Rzlbd = [cos(lbd) -sin(lbd) 0; sin(lbd) cos(lbd) 0; 0 0 1];
Ryphi = [cos(phi) 0 sin(phi);0 1 0; -sin(phi) 0 cos(phi)];
Rycor = [0 0 -1;0 1 0;1 0 0];
R = (Rzlbd*Ryphi*Rycor); % From body fixed frame to inertial frame initially aligned with principal body directions
%R2= (Rzlbd*Rzplbd*Ryphi*Rypphi*Rycor); % From body fixed frame to inertial frame initially aligned with ICRF

Aast2 = R*rdd; % Gravitational acceleration on the orbiter due to the asteroid

%Aast0 = 0;
Aast0 = -muast/norm(r)^3*r; % Gravitational acceleration on the orbiter due to 0 harmonic of the asteroid

Frad = P0*R0^2*A/norm(r-rsun)^3*(r-rsun); % Acceleration due to solar radiation

Asun = -G*Msun/norm(r-rsun)^3*(r-rsun); % Gravitation acceleration on the orbiter due to the Sun
Asunast = -G*Msun/norm(-rsun)^3*(-rsun); % Gravitational acceleration on the asteroid due to the Sun

A = (Frad/Msc + Asun-Asunast + Aast0+Aast2); % Acceleration over the orbiter

xdot=[x(4); x(5); x(6); A(1); A(2); A(3)];

end
