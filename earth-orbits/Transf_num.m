clear all;

% This code computes different types of orbit around the Earth
% Uncomment Om, w, inc, e, T and a to see to activate the different Orbits

we = 2*pi/3600/24; % Angular velocity of Earth
Re = 6400e3;            % Radius of Earth
mue= 398600.4418e9;% Earth gravitational parameter

% ------------------------------------------------ Here we define orbits

% LEO Orbit
% Om = 0/180*pi;    % Longitude of ascending node
% w  = 0/180*pi;    % Argument of periapsis
% inc= 0/180*pi;    % Inclination
% e  = 0;         % Eccentricity
% a = ((5400/2/pi)^2*mue)^(1/3); % Semi-major axis for LEO

% Transfer orbit
% Om = 0/180*pi;    % Longitude of ascending node
% w  = 0/180*pi;    % Argument of periapsis
% inc= 28/180*pi;    % Inclination
% e  = 0.7278;         % Eccentricity
% a = 0.5*((((5400/2/pi)^2*mue)^(1/3))+((mue/we^2)^(1/3))); % Semi-major axis for transfer orbit

% GEO Orbit
% Om = 0/180*pi;    % Longitude of ascending node
% w  = 0/180*pi;    % Argument of periapsis
% inc= 0/180*pi;    % Inclination
% e  = 0;         % Eccentricity
% a  = (mue/we^2)^(1/3); % Semi-major axis for GEO

% Molnyia Orbit
Om = 0;    % Longitude of ascending node
w  = 90/180*pi;    % Argument of periapsis
inc= -30/180*pi;    % Inclination
e  = 0.7;          % Eccentricity
T = 12*3600;       % Period of 12 h
a  = (mue*T^2/(2*pi)^2)^(1/3); % Semi-major axis for Molnyia

% --------------------------------------------------------------------

Rxi  = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)]; % Rotating about x for inclination

Ef = 360; % Final eccentric anomaly
xyz = zeros(3,Ef); % Vector containing cartesian coordinates
G   = zeros(3,Ef); % Ground track
for j=1:Ef
    E(j) = j/180*pi; % Our independent variable is the eccentric anomaly
    v(j) = 2*atan(sqrt((1+e)/(1-e))*tan(E(j)/2)); % Calculate true anomaly as fctn of E
    %t(j) = sqrt(a^3/mue)*(E(j)-e*sin(E(j))); % Calculate time as fctn of E
    t(j)=0; % Set time to zero for ECI frame track
    
    r = [a*(1-e^2)/(1+e*cos(v(j))); 0; 0]; % Vector pointing to the spacecraft in local coordinates
    RzOm = [cos(Om-we*t(j)) -sin(Om-we*t(j)) 0; sin(Om-we*t(j)) cos(Om-we*t(j)) 0; 0 0 1]; % Rotating about z Omega and Earth rotation
    %RzOm = [cos(Om) -sin(Om) 0; sin(Om) cos(Om) 0; 0 0 1]; % Rotating about z Omega and Earth rotation
    Rzwv = [cos(w+v(j)) -sin(w+v(j)) 0; sin(w+v(j)) cos(w+v(j)) 0; 0 0 1]; % Rotation about z of w and v
    R = RzOm*Rxi*Rzwv; % Rotation matrix
    xyz(:,j) = R'*r; % Transform to cartesian coordinates
    G(:,j) = Re/norm(xyz(:,j))*xyz(:,j); % Ground track
    lbd(j) = atan2(G(2,j),G(1,j));
end

figure(1);
hold on;
plot3([0 Re*1.1 0 0 0 0],[0 0 0 Re*1.1 0 0],[0 0 0 0 0 Re*1.1],'k');
%plot3([0 xyz(1,1)],[0 xyz(2,1)],[0 xyz(3,1)]);
plot3(xyz(1,:),xyz(2,:),xyz(3,:),'--k'); % Plot spacecraft position
plot3(G(1,:),G(2,:),G(3,:)); % Plot ground track
xlabel('x'); ylabel('y'); zlabel('z');
%camorbit(-45,25)
[xs,ys,zs] = sphere(100);
surf(Re*xs,Re*ys,Re*zs,'EdgeColor','none'); % Draw Earth
colormap copper;
alpha(0.5);
axis([-5e7 5e7 -5e7 5e7 -5e7 5e7]);

% figure(2);
% xx = linspace(0,360,Ef);
% plot(xx,lbd);




