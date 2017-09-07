% Tracking bird's flight using Kalman filter
% In this example, an observer tracks the trajectory defined by a bird.
% We assume that the bird is following a 1D accelerated motion.
% https://en.wikipedia.org/wiki/Kalman_filter contains the formulation for
% a very similar example

clear all

%% Main variables
duration = 15;  % for how long the bird flies [s]
dt = 0.1;  % observer sampling rate [s]
u = 1.5; % acceleration magnitude
X0 = [0; 0]; % initial state of the bird
BirdAccel_noise_mag = 0.5;  % process noise: the variability in how fast the bird is speeding up (stdv of acceleration: meters/sec^2)
Measurement_noise_mag = 10;   % measurement noise (stdv of location, in meters)
Ez = Measurement_noise_mag^2; % covariance of the noise (stdv^2)
Ex = BirdAccel_noise_mag^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]; % Ex converts the process noise (stdv) into covariance matrix
P = Ex; % estimate of initial bird position variance (covariance matrix)


%% System dyanmics
% An accelerated trajectory is defined by the discrete state transition
% matrix as follows:
% X = [x;v] (state)
% y = C*X (system output, or measurement)
% X(dt+1)= A*X(dt) + B*u + w (where w is the bird's acceleration noise)

A = [1 dt; 0 1] ; % state transition matrix (state prediction of the bird's flight)
B = [dt^2/2; dt]; % input control matrix
C = [1 0]; % measurement matrix (the observer only estimates position, not velocity)

% initialize output
X = X0;
array_length = floor(duration/dt) + 1;
X_loc = zeros(array_length,1); % Actual bird's flight path
X_loc_meas = zeros(array_length,1); % Bird's path as seen by the observer
vel = zeros(array_length,1); % Actual bird's velocity


%% simulate what the observer sees over time
for t_i = 1 : array_length
    
    % Generate the bird's flight
    BirdAccel_noise = BirdAccel_noise_mag * [(dt^2/2)*randn; dt*randn]; % actual acceleration, with the added noise
    X = A * X + B * u + BirdAccel_noise;
    % Generate what the observer sees
    Measurement_noise = Measurement_noise_mag * randn;
    y = C * X + Measurement_noise;
    % add new value in output arrays:
    X_loc(t_i) = X(1);
    X_loc_meas(t_i) = y;
    vel(t_i) = X(2);
    
end


%% Kalman filtering

%initialize estimation variables
X_estimate = X0;  % x_estimate of initial location estimation of where the bird is (it will be updated)
X_loc_estimate = zeros(array_length,1); %  bird position estimate
vel_estimate = zeros(array_length,1); % bird velocity estimate
P_estimate = P;
P_mag_estimate = zeros(array_length,1);
predic_state = zeros(array_length,1);
predic_var = zeros(2,2,array_length);

for t = 1:array_length
    
    % Predict next state of the bird with the last state and predicted motion.
    X_estimate = A * X_estimate + B * u;
    predic_state(t) = X_estimate(1);
    %predict next covariance
    P = A * P * A' + Ex;
    predic_var(:,:,t) = P;
    % predicted observer's measurement covariance
    % Kalman Gain
    K = P*C'*inv(C*P*C'+Ez);
    % Update the state estimate.
    X_estimate = X_estimate + K * (X_loc_meas(t) - C * X_estimate);
    % update covariance estimation.
    P =  (eye(2)-K*C)*P;
    %Store for plotting
    X_loc_estimate(t) = X_estimate(1);
    vel_estimate(t) = X_estimate(2);
    P_mag_estimate(t) = P(1);
    
end


% plot location
figure(1);
clf;
hold on;
x_axis = 0:dt:duration;
title('Position');
plot(x_axis, X_loc, '-r.');
plot(x_axis, X_loc_meas, '-k.');
plot(x_axis, X_loc_estimate,'-g.');
legend('Actual','Observed','Estimated');
axis([0 duration min(min(X_loc_meas),min(X_loc))-10 max(max(X_loc_meas),max(X_loc))+10]);

% plot velocity
figure(2);
clf;
hold on;
title('Velocity');
plot(x_axis, vel, '-r.');
plot(x_axis, vel_estimate, '-g.');
legend('Actual','Estimated');
axis([0 duration min(vel)-1 max(vel)+1]);

