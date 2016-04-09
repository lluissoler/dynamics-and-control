clear all;

R0 = 149597870000; % 1 AU [m]
Rs = 2.497*R0; % Asteroid orbit radius [m]
G  = 6.673e-11; % Gravitational constant [m/kg/s2]
Msc= 800; % Mass of the orbiter [kg]
Msun = 1.9891e30; % Mass of the Sun [kg]
%Msun = 0; % Mass of the Sun [kg]
muast= 14.043; % Gravitational parameter of the asteroid
i = 2.277/180*pi; % Asteroid inclination orbit angle [rad]
plbd = 202/180*pi; % Asteroid's pole longitude [rad]
pphi = -45/180*pi; % Asteroid's pole latitude [rad]
Ts = sqrt(Rs^3/G/Msun)*2*pi; % Period of the asteroid orbit [s]
W = 1/(Ts/2/pi); % Angular velocity of the asteroid's orbit [rad/s]
Ta = 6.02978*3600; % Rotational period of the asteroid [s]
w = 1/(Ta/2/pi); % Angular velocity of the asteroid [rad/s]
P0 = 4.56e-6; % Solar pressure at 1 AU [N/m2]
%P0 = 0; % Solar pressure at 1 AU [N/m2]
A  = 40; % Orbiter surface [m2]
C20 = -0.0142503; % C20 spherical harmonic
C22 = 0.0229886; % C22 spherical harmonic
Rref = 396.1; % Reference radius of the asteroid

T = 3600*24*30; % 1 month [s]
ic = 3; % Sets the initial condition to be taken
if ic==1
    x0 = [551.36; 0; 0; 0; 0; 0]; % Initial condition 1, body fixed frame
elseif ic==2
    x0 = [0; 551.36; 0; -0.1596; 0; 0]; % Initial condition 2, body aligned frame
elseif ic==3
    x0 = [0; 0; 551.36; 0; 0.1596; 0]; % Initial condition 3, ICRF aligned frame
end

% Solving dynamics
RelTol = 1e-6;
AbsTol = 1e-6;
options = odeset('RelTol',RelTol,'AbsTol',[AbsTol AbsTol AbsTol AbsTol AbsTol AbsTol]);
[t,x]=ode45(@(t,x) f_orbit(t,x,R0,Rs,G,Msc,Msun,muast,i,W,P0,A,w,plbd,pphi,C20,C22),[0 T],x0,options);
last = length(t);

r = [x(:,1) x(:,2) x(:,3)]; % Position of the spacecraft

% for i=1:last % Rotate the spacecraft when body fixed frame
%     Rzlbdt = [cos(-w*t(i)) -sin(-w*t(i)) 0; sin(-w*t(i)) cos(-w*t(i)) 0; 0 0 1];
%     r(i,:)= (Rzlbdt*r(i,:)')';
% end


modr = zeros(size(t));
for i=1:last
    modr(i) = sqrt(r(i,1)^2+r(i,2)^2+r(i,1)^2);
end

rsun = [-Rs*cos(W*t) -Rs*sin(W*t)*cos(i) Rs*sin(W*t)*sin(i)]; % Position of the Sun
% Rzplbd= [cos(plbd) -sin(plbd) 0; sin(plbd) cos(plbd) 0; 0 0 1];
% Rypphi= [cos(pphi) 0 sin(pphi);0 1 0; -sin(pphi) 0 cos(pphi)];
% rsun = (Rzplbd*Rypphi*rsun')'; % Rotating position of the Sun if our inertial frame is the one in P1

phi = pi/2-asin(r(:,3)./modr(:));
lbd = atan2(r(:,2),r(:,1));

for i=1:last
    Asun(i,:) = -G*Msun/norm(r(i,:)-rsun(i,:))^3*(r(i,:)-rsun(i,:)); % Gravitation force on the orbiter due to the Sun
    Asunast(i,:) = -G*Msun/norm(-rsun(i,:))^3*(-rsun(i,:));
    Asunrel(i,:) = Asun(i,:)-Asunast(i,:);
        
    Fsun(i,:) = Asunrel(i,:)*Msc;
    Frad(i,:) = P0*R0^2*A/norm(r(i,:)-rsun(i,:))^3*(r(i,:)-rsun(i,:)); % Force due to solar radiation
    Fast(i,:) = -muast*Msc/norm(r(i,:))^3*r(i,:);
    
    modFrad(i) = sqrt(Frad(i,1)^2+Frad(i,2)^2+Frad(i,1)^2);
    modFsun(i) = sqrt(Fsun(i,1)^2+Fsun(i,2)^2+Fsun(i,1)^2);
    modFast(i) = sqrt(Fast(i,1)^2+Fast(i,2)^2+Fast(i,1)^2); 
end
% 
figure(1);
hold on;
axis([-800,800,-800,800,-800,800]);
xlabel('x');
ylabel('y');
zlabel('z');
plot(0,0,'x');
text(50,50,'Golevka');
plot3(r(:,1),r(:,2),r(:,3));
[xs,ys,zs] = sphere(100);
surf(Rref*xs,Rref*ys,Rref*zs,'EdgeColor','none'); % Draw asteroid
colormap copper;
alpha(0.5);
% plot(rsun(:,1),rsun(:,2));
% 
% 
% figure(2);
% hold on;
% axis([-1000,1000,-1000,1000]);
% for i=1:25
%     k = i*15;
%     cla;
%     plot([r(k,1) r(k,1)+Frad(k,1)],[r(k,2) r(k,2)+Frad(k,2)],'b',[r(k,1) r(k,1)+Fsun(k,1)],[r(k,2) r(k,2)+Fsun(k,2)],'r',[r(k,1) r(k,1)+Fast(k,1)],[r(k,2) r(k,2)+Fast(k,2)],'g',[0 0],[0 0],'x',[r(k,1) r(k,1)],[r(k,2) r(k,2)],'x');
%     plot([r(k,1) r(k,1)+Frad(k,1)],[r(k,2) r(k,2)+Frad(k,2)],[r(k,1) r(k,1)+Fsun(k,1)],[r(k,2) r(k,2)+Fsun(k,2)],[r(k,1) r(k,1)+Fast(k,1)],[r(k,2) r(k,2)+Fast(k,2)],[0 0],[0 0],'x',[rsun(k,1) rsun(k,1)],[rsun(k,2) rsun(k,2)],'o');
%     %legend('Frad','Fsun','Fast');
%     plot([r(k,1) rsun(k,1)],[r(k,2) rsun(k,2)],'c:',[0 r(k,1)],[0 r(k,2)],'c:');
%     plot([r(k,1) r(k,1)+r(k,4)*10],[r(k,2) r(k,2)+r(k,5)*10],'y-.');
%     text(-900,-900,['T=' num2str(t(i))]);
%     pause(0.5);
% end

% figure(3);
% plot(t,r(:,4),t,r(:,5),t,r(:,6));
% % %plot(t,r(:,5));
% % %plot(t,modr,t,modFast'.*modr.^3);


figure(3);
hold on;
% [xs,ys,zs] = sphere(100);
% surf(Rref*xs,Rref*ys,Rref*zs,'EdgeColor','none'); % Draw asteroid
% colormap copper;
% alpha(0.5);
axis([-800,800,-800,800,-800,800]);
xlabel('x');
ylabel('y');
zlabel('z');
plot(0,0,'x');
text(50,50,'Golevka');
[xs,ys,zs] = sphere(100);
surf(Rref*xs,Rref*ys,Rref*zs,'EdgeColor','none'); % Draw asteroid
colormap copper;
alpha(0.5);
plot3([0 rsun(1,1)],[0 rsun(1,2)],[0 rsun(1,3)],'k:'); % Initial Sun position
%plot3([0 rsun(last,1)],[0 rsun(last,2)],[0 rsun(last,3)],'k:'); % Final Sun position
for i=1:200
    ri(i,:) = r(i,:);
    rf(i,:) = r(last-i,:);
end
plot3(ri(:,1),ri(:,2),ri(:,3),'b-');
plot3(rf(:,1),rf(:,2),rf(:,3),'r-.');
%legend('','','','Initial orbit (t=0)','Final orbit (t=1 month)');

%
%  figure(3);
%  hold on;
%  plot(t,Frad);
 
 
 
 