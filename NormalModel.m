% Initialization
clear all;

deltat = 0.001;
xd = 1;    % Typical relaxed fiber length during diastole
xs = 0.5;  % Assumed contraction fiber length
bd = 1;    % Baseline parameter for diastole
bs = 1.5;  % Adjusted baseline parameter for systole, assuming it should be higher
T = 1;     % Constant proportional to tension
e = 0.001; % Small constant for timescale
N = 20000; % Number of iterations to cover more heartbeats

% Parameters for automatic pulsing based on 'b' value
b_threshold = -0.9;  % b value threshold for triggering a pulse
pulseWidth = 100;    % Width of the pulse in iteration steps, simulate the pulse duration
pulseFrequency = 1000; % Frequency of pulses in terms of iteration steps
AV_decay_rate = 0.005; % Rate of potential decrease when no signal received
inPulse = false;
atrialSignalCounter = 0; % Initialize counter for atrial signals

% Initial values
x = zeros(1, N+1);
b = zeros(1, N+1);
u_values = zeros(1, N+1); % To store the pulse state for visualization
x(1) = 0.8;
b(1) = -0.9;
u = 0;

% Simulation with automatic pulsing based on 'b' value
for i = 1:N
    t = i * deltat; % Current time step
    
    % Check if it's time to start a new pulse based on frequency
    if mod(i, pulseFrequency) == 0
        atrialSignalCounter = atrialSignalCounter + 1; % Increment atrial signal counter
        % Check if 'b' falls below the threshold to trigger a pulse
        if b(i) < b_threshold && ~inPulse
            inPulse = true;
            pulseCounter = pulseWidth; % Reset pulse counter
            if mod(atrialSignalCounter, 2) == 0
                pulseAmplitude = 0.5; % Reduce the amplitude for every other pulse
            else
                pulseAmplitude = 1; % Normal pulse
            end
        end
    end
    
    % Generate pulse if in pulse period
    if inPulse
        u = pulseAmplitude; % Use the modified amplitude
        pulseCounter = pulseCounter - 1;
        if pulseCounter <= 0
            inPulse = false; % End of pulse
        end
    else
        u = 0;
    end
    
    u_values(i) = u; % Store the current value of u for visualization
    
    % Update b based on the current state and presence of a pulse
    b(i+1) = b(i) - AV_decay_rate * deltat + deltat * ((x(i) - xd) + 10 * (xd - xs) * u);
    
    % Update x with added non-linearity for realism
    x(i+1) = x(i) + deltat * ((-1/e) * (x(i)^3 - x(i)^2 + T * x(i) + b(i)));
end

% Plot results
t_start = 0;  % Start time
t_end = 20;   % End time
idx_start = round(t_start / deltat) + 1;
idx_end = round(t_end / deltat) + 1;

figure(1)
subplot(2,1,1);
plot((idx_start:idx_end) * deltat, x(idx_start:idx_end), 'LineWidth', 2);
title('Fiber Length x(t)');
xlabel('Time');
ylabel('x');

subplot(2,1,2);
plot((idx_start:idx_end) * deltat, b(idx_start:idx_end), 'LineWidth', 2);
title('Electrical Control Variable b(t)');
xlabel('Time');
ylabel('b');


%subplot(3,1,3); % Plotting the pulse signal
%plot((idx_start:idx_end) * deltat, u_values(idx_start:idx_end), 'LineWidth', 2);
%title('Pulse Signal u(t)');
%xlabel('Time');
%ylabel('u');
