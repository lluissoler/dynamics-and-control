function xdot=f_orbit(x,t,Msat,mue,mum,P0,A,Reb,Rmb,e,phi,Sol_thr)

r = [x(1); x(2); x(3)]; % Position vector of the spacecraft
rd= [x(4); x(5); x(6)]; % Velocity vector of the spacecraft

Re_sc = r-Reb; % Vector from center of Earth to satellite
% if (Re_sc)'*rd >= 0 % Calculates the true anomaly of the orbit. Note that this will only work for eliptical orbits
%     nu = acos([-e*cos(phi) 0 e*sin(phi)]*Re_sc/e/norm(Re_sc));
% else
%     nu = 2*pi - acos([-e*cos(phi) 0 e*sin(phi)]*Re_sc/e/norm(Re_sc));
% end

R0 = 149598000000; % 1 AU
rho = 0.9;
wSun = 2*pi/365.25/24/60/60; % Angular velocity of the Sun around the Earth
wMoon= 2*pi/27.32/24/60/60; % Angular velocity of the Earth-Moon system
%wMoon=0;
iem = 5.145/180*pi; % Inclination of Eart-Moon system
Rsun = [-R0*cos((wSun+wMoon)*t)*cos(iem); -R0*sin((wSun+wMoon)*t); R0*cos((wSun+wMoon)*t)*sin(iem)]+Reb; % Position of the Sun wrt our frame

Ae = -mue/norm(Re_sc)^3*(Re_sc); % Gravitational acceleration of the satellite due to Earth's gravitational field
%Am = -mum/norm(r-Rmb)^3*(r-Rmb); % Gravitational acceleration of the satellite due to Moon's gravitational field
Am = 0; % If not considering Moon's gravity field

SLim = 15/180*pi; % 180+-SLim is the range where solar pressure will be used as thrust
% if Sol_thr == 1 && nu > pi-SLim && nu < pi+SLim
%     Arad = 5.89e-5*rd/norm(rd); % Acceleration due to solar radiation [m/s2], opposite to the velocity (tangential)
%     %Arad = 0;
% else
%     Arad = 0;
% end
if Sol_thr == 1 && (r-Rsun)'*rd > 0
    % Arad = 5.89e-5*rd/norm(rd)/sqrt(2); % Old value
    % Arad = 1e-5*rd/norm(rd); % Using the mean value at all times
    % c_alpha = (r-Rsun)'*rd/norm(r-Rsun)/norm(rd);
    % Arad = P0*R0^2*A/norm(r-Rsun)^2/Msat*c_alpha*rho*rd/norm(rd); % Solar pressure pointing to the direction of v (rd)
    Arad = P0*R0^2*A/norm(r-Rsun)^3/Msat*rho*(r-Rsun); % Solar pressure always on the direction from the Sun (rd)
else
    Arad = 0;
    % c_alpha = 0;
end

Acc = (Arad + Ae + Am); % Total Acceleration over the orbiter
t*3.1710e-008
xdot=[x(4); x(5); x(6); Acc(1); Acc(2); Acc(3)];

end
