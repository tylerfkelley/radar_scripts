% Parameters
fs = 3e9; % Sampling frequency (Hz)
f0 = 1e9; % Start frequency of chirp (Hz)
bw = 200e6; % Chirp bandwidth (Hz)
T = 1e-3; % Chirp duration (s)

% Target parameters
target1Distance = 10; % Target 1 at 10 meters
target1Speed = 0;    % Target 1 moving at 50 m/s
target2Distance = 10; % Target 2 at 30 meters
target2Speed = -50;   % Target 2 moving at -20 m/s (approaching)

c = 3e8; % Speed of light (m/s)

% Control flag for plotting waveforms
plotWaveforms = false;

% Time vector
t = (0:1/fs:T-1/fs);

% Generate transmitted chirp signal
chirpRate = bw / T;
txSignal = exp(1j * 2 * pi * (f0 * t + (chirpRate/2) * t.^2));

%% --- Compute Reflected Signal for Target 1 ---
timeDelay1 = 2 * target1Distance / c;
instantaneousFreq1 = f0 + chirpRate * t;
dopplerShift1 = 2 * instantaneousFreq1 * target1Speed / c;

delayedSignal1 = interp1(t, txSignal, t - timeDelay1, 'linear', 0);
reflectedSignal1 = delayedSignal1 .* exp(1j * 2 * pi * cumsum(dopplerShift1) / fs);

% IQ Mixer Processing for Target 1
A1 = conj(reflectedSignal1) .* txSignal;
I1 = real(A1);
Q1 = imag(A1);

%% --- Compute Reflected Signal for Target 2 ---
timeDelay2 = 2 * target2Distance / c;
instantaneousFreq2 = f0 + chirpRate * t;
dopplerShift2 = 2 * instantaneousFreq2 * target2Speed / c;

delayedSignal2 = interp1(t, txSignal, t - timeDelay2, 'linear', 0);
reflectedSignal2 = delayedSignal2 .* exp(1j * 2 * pi * cumsum(dopplerShift2) / fs);

% IQ Mixer Processing for Target 2
A2 = conj(reflectedSignal2) .* txSignal;
I2 = real(A2);
Q2 = imag(A2);

%% --- Optional Plot of Transmit and Received Waveforms ---
if plotWaveforms
    figure;
    subplot(2,1,1);
    plot(t * 1e3, real(txSignal));
    title('Transmitted Signal');
    xlabel('Time (ms)');
    ylabel('Amplitude');
    grid on;
    
    subplot(2,1,2);
    plot(t * 1e3, real(reflectedSignal1), t * 1e3, real(reflectedSignal2));
    title('Received Signals');
    xlabel('Time (ms)');
    ylabel('Amplitude');
    legend('Target 1', 'Target 2');
    grid on;
end

%% --- Plot Results ---
if plotWaveforms
    figure
else
    clf;
end

% 1. Real and Imaginary parts of mixed signal 1
subplot(3,1,1);
plot(t * 1e3, I1, t * 1e3, Q1);
title('Target 1 - IQ Mixer Output');
xlabel('Time (ms)');
ylabel('Amplitude');
legend('In-phase (I)', 'Quadrature (Q)');
grid on;

% 2. Real and Imaginary parts of mixed signal 2
subplot(3,1,2);
plot(t * 1e3, I2, t * 1e3, Q2);
title('Target 2 - IQ Mixer Output');
xlabel('Time (ms)');
ylabel('Amplitude');
legend('In-phase (I)', 'Quadrature (Q)');
grid on;

% 3. Phase of both signals
subplot(3,1,3);
plot(t * 1e3, (angle(A1)), t * 1e3, (angle(A2)));
title('Phase of Target 1 and Target 2');
xlabel('Time (ms)');
ylabel('Phase (radians)');
legend('Target 1', 'Target 2');
grid on;