clear all;

% ------
% Inputs
% ------

deltaV1 = 120;
deltaV2 = 7.6;

Me = 5.9736e24; % Mass of the Earth
Mm = 7.3477e22; % Mass of the Moon
G  = 6.673e-11; % Gravitational constant
%iem= 5.145/180*pi; % Inclination of Eart-Moon system
%Eat= 23.44/180*pi; % Earth axial tilt
iem = 0; Eat = 0;

Rem = 384399e3; % Semi-major axis of Moon's orbit around the Earth
mue = G*Me; % Gravitational parameter of the Earth
mum = G*Mm; % Gravitational parameter of the Moon

Rb = (Rem*Mm)/(Me+Mm); % Distance of the Earth-Moon barycenter from the center of the Earth
Reb= [-Rb; 0; 0];
Rmb= [Rem-Rb; 0; 0];

R0 = 149598000000; % 1 AU
rho = 0.9; % Reflectivity
wSun = 2*pi/365.25/24/60/60; % Angular velocity of the Sun around the Earth
wMoon= 2*pi/27.32/24/60/60; % Angular velocity of the Earth-Moon system

Re = 6371e3; % Radius of the Earth
Rm = 1737.1e3; % Radius of the Moon
P0 = 4.56e-6; % Solar pressure at 1 AU [N/m2]
Msat = 9.5585; % Mass of satellite
A = 78; % Effective area of the solar sail
rho = 0.9; % Reflectivity

% ------------------
% Orbit calculations
% ------------------

% GTO orbit
Rper = 6678e3; % GTO Periapse
Rapo = 42164e3;% GTO Apoapse = GEO radius
a = 0.5*(Rapo+Rper); % GTO semi-major axis
e = 0.5*(-Rper/a + sqrt((Rper/a)^2-4*(Rper/a-1))); % GTO eccentricity
p = a*(1-e^2);
v0= sqrt(mue/p)*sqrt(1+2*e+e^2); % Velocity at periapse of GTO
va= sqrt(mue/p)*sqrt(1-2*e+e^2); % Velocity at apoapse of GTO

% GEO orbit: Find necessary deltaVGEO at apoapse of GTO to switch to GEO
aGEO = Rapo; pGEO = Rapo; % GEO
vGEO = sqrt(mue/pGEO); % Velocity at periapse of GEO
deltaVGEO = vGEO - va;

% -----------
% Integration
% -----------
RelTol = 1e-7;
AbsTol = 1e-7;
options = odeset('RelTol',RelTol,'AbsTol',[AbsTol AbsTol AbsTol AbsTol AbsTol AbsTol]);
ini_orb = 'GTO';

% Initial GTO/GEO orbit
% ---------------------
if ini_orb == 'GTO'
    T1 = 2*pi*sqrt(a^3/mue); % Final time 1 (one complete period of GTO orbit)
    x0 = [-Rper*cos(iem+Eat)-Rb; 0; Rper*sin(iem+Eat); 0; -v0; 0]; % Initial conditions of the satellite
elseif ini_orb == 'GEO'
    T1 = 2*pi*sqrt(aGEO^3/mue); % Final time 2 (period of GEO orbit)
    x0 = [-Rapo*cos(iem+Eat)-Rb; 0; Rapo*sin(iem+Eat); 0; -vGEO; 0]; % Initial conditions of the satellite
end
[t1,x1]=ode45(@(t,x) f_orbit(x,t,Msat,mue,mum,P0,A,Reb,Rmb,e,iem+Eat,0),[0 T1],x0,options);
last1 = length(t1);

for i=1:last1-1 % We do not include the last row of this matrix, since this state is substituted by first row of next orbit
    r(i,:) = [x1(i,1) x1(i,2) x1(i,3)]; % Position of the spacecraft
    mod_r(i) = norm(r(i,:)); % Distance from the origin
    v(i,:) = [x1(i,4) x1(i,5) x1(i,6)]; % Velocity of the spacecraft
%     if (r(i,:)-Reb')*v(i,:)' >= 0
%         nu(i) = acos([-e 0 0]*(r(i,:)-Reb')'/e/norm(r(i,:)-Reb'));
%     else
%         nu(i) = 2*pi - acos([-e 0 0]*(r(i,:)-Reb')'/e/norm(r(i,:)-Reb'));
%     end
end


% First deltaV at periapse
% ------------------------

% Transfer orbit (first deltaV at periapse)
T2 = 60*60*24*1; % Final time 2 (half complete period of transfer orbit)
x0 = [x1(last1,1); x1(last1,2); x1(last1,3); x1(last1,4); x1(last1,5)-deltaV1; x1(last1,6)]; % Initial conditions of the satellite
[t2,x2]=ode45(@(t,x) f_orbit(x,t,Msat,mue,mum,P0,A,Reb,Rmb,e,iem+Eat,0),[0 T2],x0,options);
last2 = length(t2);

for i=1:last2-1
    r(last1-1+i,:) = [x2(i,1) x2(i,2) x2(i,3)]; % Position of the spacecraft
    v(last1-1+i,:) = [x2(i,4) x2(i,5) x2(i,6)]; % Velocity of the spacecraft
    
    x21(i,:) = [x2(i,1) x2(i,2) x2(i,3)]-Reb';
    mod_r2(i) = norm(x21(i,:)); % Distance from the origin
    
    mod_r(last1-1+i) = norm(r(last1-1+i,:)); % Distance from the origin
%     if (r(last1-1+i,:)-Reb')*v(last1-1+i,:)' >= 0
%         nu(last1-1+i) = acos([-e 0 0]*(r(last1-1+i,:)-Reb')'/e/norm(r(last1-1+i,:)-Reb'));
%     else
%         nu(last1-1+i) = 2*pi - acos([-e 0 0]*(r(last1-1+i,:)-Reb')'/e/norm(r(last1-1+i,:)-Reb'));
%     end
end

[max_dist pos_max_dist] = max(mod_r2);
x22 = zeros(pos_max_dist,3);
for i=1:pos_max_dist
    x22(i,:) = [x2(i,1) x2(i,2) x2(i,3)];
end

% Final orbit (second deltaV at apoapse)
T3 = 60*60*24*50; % Final time 2
x0 = [x22(pos_max_dist,1); x22(pos_max_dist,2); x22(pos_max_dist,3); x2(pos_max_dist,4); x2(pos_max_dist,5)+deltaV2; x2(pos_max_dist,6)]; % Initial conditions of the satellite
[t3,x3]=ode45(@(t,x) f_orbit(x,t,Msat,mue,mum,P0,A,Reb,Rmb,0,iem+Eat,1),[0 T3],x0,options);
last3 = length(t3);

for i=1:last3-1
    r(pos_max_dist-2+i,:) = [x3(i,1) x3(i,2) x3(i,3)]; % Position of the spacecraft
    v(pos_max_dist-2+i,:) = [x3(i,4) x3(i,5) x3(i,6)]; % Velocity of the spacecraft
%     if [-1 0 0]*v(pos_max_dist-2+i,:)' <= 0
%         nu(pos_max_dist-2+i) = acos([-1 0 0]*(r(pos_max_dist-2+i,:)-Reb')'/norm(r(pos_max_dist-2+i,:)-Reb'));
%     else
%         nu(pos_max_dist-2+i) = 2*pi - acos([-1 0 0]*(r(pos_max_dist-2+i,:)-Reb')'/norm(r(pos_max_dist-2+i,:)-Reb'));
%     end
end

Rsun = [-R0*cos((wSun+wMoon)*t3)*cos(iem)+Reb(1) -R0*sin((wSun+wMoon)*t3)+Reb(2) R0*cos((wSun+wMoon)*t3)*sin(iem)+Reb(3)]; % Position of the Sun wrt our frame
r3 = [x3(:,1) x3(:,2) x3(:,3)];
v3 = [x3(:,4) x3(:,5) x3(:,6)];
c_alpha = zeros(last3,1);
Arad = zeros(last3,3);
mod_Arad = zeros(1,last3);
for i=1:last3
    E(i) = 0.5*norm(v3(i,:))^2 - mue/norm(r3(i,:)-Reb');
    if (r3(i,:)-Rsun(i,:))*v3(i,:)' > 0
        % Arad = 5.89e-5*rd/norm(rd)/sqrt(2);
        c_alpha(i) = (r3(i,:)-Rsun(i,:))*v3(i,:)'/norm(r3(i,:)-Rsun(i,:))/norm(v3(i,:));
        
        Arad1(i,:) = P0*R0^2*A/norm(r3(i,:)-Rsun(i,:))^2/Msat*c_alpha(i)*rho*v3(i,:)/norm(v3(i,:)); % Solar pressure pointing to the direction of v (rd)
        mod_Arad1(i) = norm(Arad1(i,:));
        
        Arad2(i,:) = P0*R0^2*A/norm(r3(i,:)-Rsun(i,:))^3/Msat*rho*(r(i,:)-Rsun(i,:)); % Solar pressure pointing to the direction of v (rd)
        mod_Arad2(i) = norm(Arad2(i,:));
        
        %Arad = 1e-5*rd/norm(rd);
    else
        Arad1(i,:) = [0 0 0];
        mod_Arad1(i) = 0;
        Arad2(i,:) = [0 0 0];
        mod_Arad2(i) = 0;
    end
    i/last3
end

% figure(2);
% plot(t3,mod_Arad1,t3,mod_Arad2);
% legend('Parallel to v','Parallel to Rsun');

plot(t3,E);

% --------
% Plotting
% --------
[xs,ys,zs] = sphere(25); % Draw the Earth

figure(1);
hold on;
xlabel('x'); ylabel('y'); zlabel('z');
plot3(0,0,0,'x');
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]);
surf((Re*xs-Rb)/Rem,Re*ys/Rem,Re*zs/Rem,'EdgeColor','none'); % Draw the Earth
colormap winter; % Draw the Earth
surf((Rm*xs+(Rem-Rb))/Rem,Rm*ys/Rem,Rm*zs/Rem,'EdgeColor','none'); % Draw the Moon
colormap gray; % Draw the Moon
alpha(0.5); % Draw the Moon

plot3(x1(:,1)/Rem,x1(:,2)/Rem,x1(:,3)/Rem,'b'); % Plot the trajectory of the satellite (GTO orbit)
plot3(x22(:,1)/Rem,x22(:,2)/Rem,x22(:,3)/Rem,'r'); % Plot the trajectory of the satellite (transfer orbit)
plot3(x3(:,1)/Rem,x3(:,2)/Rem,x3(:,3)/Rem,'g'); % Plot the trajectory of the satellite (transfer orbit)
hold on;

%plot3(r(:,1)/Rem,r(:,2)/Rem,r(:,3)/Rem);









