% This is an FMCW simulation and also to practice OOP
clear; close all;

addpath("utils_functions_kelley\");

FONTSIZE = 20;
gca_linewidth = 2.5;
set(0, 'defaultfigurecolor', [1 1 1]); % White background for plots
set(0, 'defaultAxesFontSize', FONTSIZE);
set(0, 'defaultlinelinewidth', 2);
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');

c = 3e8; % Speed of light

% Define start and stop frequencies of the waveform
f_start = 0.4e9; % 0.4 GHz
f_stop = 4e9; % 4 GHz
T_chirp = 1e-6; % Chirp duration (1 microsecond)
Fs = 256e9; % Sample frequency
NumPoints = round(Fs * T_chirp); % Number of samples per chirp

% Number of chirps to process
N_chirps = 1; % Set this to control how many chirps to process

% Generate extended time vector for multiple chirps
T_total = N_chirps * T_chirp; % Total time duration
NumTotalPoints = N_chirps * NumPoints; % Total number of points

t = (0:NumTotalPoints - 1) / Fs; % Extended time vector

% Compute bandwidth and chirp slope
B = f_stop - f_start;
slope = B / T_chirp;
fft_resolution = 1 / T_chirp;


% Generate chirp signal for N_chirps
phase_chirp = 2 * pi * (f_start * t + (B / (2 * T_chirp)) * mod(t, T_chirp).^2);
s_t = cos(phase_chirp); % Transmitted signal
s_t = hilbert(s_t);


% Define target distances
d1 = 0.03; % Target 1 at 2.2 meters
d2 = 0.2; % Target 2 at 2.30 meters
d3 = 0.35; % Target 2 at 2.30 meters

disp(['Distance Between Targets: ', num2str(abs(d2 - d1)), ' meters']);

% Compute round-trip delays
tau1 = 2 * d1 / c;
tau2 = 2 * d2 / c;
tau3 = 2 * d3 / c;

% Compute sample shifts
N_shift1 = round(tau1 * Fs);
N_shift2 = round(tau2 * Fs);
N_shift3 = round(tau3 * Fs);

% Generate received signals by shifting the transmitted chirp
s_received1 = circshift(s_t, N_shift1);
s_received2 = circshift(s_t, N_shift2);
s_received3 = circshift(s_t, N_shift3);

%We now do the hilbert transform of the received signal to get the real
%and the imaginary part of the signal.

% Superimpose the received signals
s_received = s_received1 + 0.5*s_received2 + 0.25*s_received3;

beat_signal = conj(s_received) .* s_t;







% Plot signals
figure;
subplot(3, 1, 1)
hold on;
plot(t * 1e6, real(s_t), 'b', 'LineWidth', 2);
plot(t * 1e6, imag(s_t), 'r', 'LineWidth', 2);
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Transmitted\: and\: Received\: Signals}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 0.1]) % Adjust for multiple chirps
legend ("I" , "Q")

subplot(3, 1, 2)
plot(t * 1e6, real(s_received), 'b', 'LineWidth', 2);
hold on;
plot(t * 1e6, imag(s_received), 'r', 'LineWidth', 2);
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Transmitted\: and\: Received\: Signals}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
legend ("I" , "Q")
xlim([0 0.1]) % Adjust for multiple chirps

subplot(3, 1, 3)
plot(t * 1e6, real(beat_signal), 'b', 'LineWidth', 2);
hold on;
plot(t * 1e6, imag(beat_signal), 'r', 'LineWidth', 2);
%plot(t * 1e6, beat_signal_Q + beat_signal_I, 'r', 'LineWidth', 2);
title("$$\bf{Mixed\: Signal}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 0.1]) 
legend ("I" , "Q")

% Compute FFT of the beat signal
N_fft = 2^nextpow2(NumTotalPoints); % Use total points for FFT
f_axis = linspace(0, Fs / 2, N_fft / 2); % Frequency axis

% Perform FFT and normalize

S_beat_FFT = abs(fft(real(beat_signal), N_fft)); % Compute FFT magnitude
S_beat_FFT = S_beat_FFT(1:N_fft / 2); % Keep positive frequencies
S_beat_FFT = S_beat_FFT ./ max(S_beat_FFT); % Normalize

% Convert beat frequency to range
range_axis = beat2range(f_axis', slope);


% Plot Range Spectrum
figure;
plot(range_axis, S_beat_FFT, 'b', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Normalized\: Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Range\: Spectrum\: Old \: Way}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 5]);

% Plot Range Spectrum
figure;
stem(range_axis, S_beat_FFT, 'b', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Normalized\: Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Range\: Spectrum\: Actual \: Points}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 5]);

%%

S_beat_FFT_interp = interpft(S_beat_FFT,length(S_beat_FFT).*10);

% Define the new interpolated frequency axis
f_axis_interp = linspace(0, Fs / 2, length(S_beat_FFT_interp));

% Interpolate the range axis by the same factor (10x)
range_axis_interp = beat2range(f_axis_interp', slope);

% Plot Range Spectrum
figure;
plot(range_axis_interp,S_beat_FFT_interp, 'b', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Normalized\: Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Range\: Spectrum\: Up \: FFT \: Interpolation}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 5]);

figure;
hold on;
plot(range_axis_interp,S_beat_FFT_interp, 'r', 'LineWidth', 2);
stem(range_axis,S_beat_FFT, 'b', 'LineWidth', 2);
%plot(range_axis,S_beat_FFT, 'k', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Normalized\: Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Range\: Spectrum\: Up \: FFT \: Interpolation}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
legend("Interpolated FFT" , "Actual FFT Points", "Previous Fit for FFT")
xlim([0.5 1.5]);
ylim([0 1]);


%% Now we apply windowing 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform FFT and normalize

window_blackman = blackmanharris(length(real(beat_signal)));


S_beat_FFT = abs(fft(real(beat_signal) .* window_blackman', N_fft)); % Compute FFT magnitude
S_beat_FFT = S_beat_FFT(1:N_fft / 2); % Keep positive frequencies
S_beat_FFT = S_beat_FFT ./ max(S_beat_FFT); % Normalize

% Convert beat frequency to range
range_axis = beat2range(f_axis', slope);


% Plot Range Spectrum
figure;
plot(range_axis, S_beat_FFT, 'b', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Normalized\: Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Range\: Spectrum\: Old \: Way}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 5]);

% Plot Range Spectrum
figure;
stem(range_axis, S_beat_FFT, 'b', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Normalized\: Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Range\: Spectrum\: Actual \: Points :\ (WINDOWED)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 5]);

%%

S_beat_FFT_interp = interpft(S_beat_FFT,length(S_beat_FFT).*5);

% Define the new interpolated frequency axis
f_axis_interp = linspace(0, Fs / 2, length(S_beat_FFT_interp));

% Interpolate the range axis by the same factor (10x)
range_axis_interp = beat2range(f_axis_interp', slope);

% Plot Range Spectrum
figure;
plot(range_axis_interp,S_beat_FFT_interp, 'b', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Normalized\: Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Windowed\ Range\: Spectrum\: Up \: Interpolated}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 5]);

%% This is to see the zero padding method...

beat_signal_padded  = Kelley.zpad(beat_signal,10)';
NumTotalPoints_padded = length(beat_signal_padded);
% Compute FFT of the beat signal
N_fft_padded = 2^nextpow2(NumTotalPoints_padded); % Use total points for FFT
f_axis_padded = linspace(0, Fs / 2, N_fft_padded / 2); % Frequency axis


S_beat_FFT = abs(fft(real(beat_signal_padded), N_fft_padded)); % Compute FFT magnitude
S_beat_FFT = S_beat_FFT(1:N_fft_padded / 2); % Keep positive frequencies
S_beat_FFT = S_beat_FFT ./ max(S_beat_FFT); % Normalize

% Convert beat frequency to range
range_axis = beat2range(f_axis_padded', slope);


% Plot Range Spectrum
figure;
plot(range_axis, S_beat_FFT, 'k', 'LineWidth', 3);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Normalized\: Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Range\: Spectrum\: Old \: Way}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 1]);
