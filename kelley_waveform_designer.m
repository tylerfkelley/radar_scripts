% This is an FMCW simulation and also to practice OOP
clear; close all;

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
N_chirps = 5; % Set this to control how many chirps to process

% Generate extended time vector for multiple chirps
T_total = N_chirps * T_chirp; % Total time duration
NumTotalPoints = N_chirps * NumPoints; % Total number of points

t = (0:NumTotalPoints - 1) / Fs; % Extended time vector

% Compute bandwidth and chirp slope
B = f_stop - f_start;
slope = B / T_chirp;
fft_resolution = 1 / T_chirp;

disp(['Range Resolution: ', num2str(c / (2 * B)), ' meters']);
disp(['FFT Resolution: ', num2str((1 / T_chirp) / 1e6), ' MHz']);

diff_range = beat2range(2 * fft_resolution, slope) - beat2range(fft_resolution, slope);
disp(['Alternate Range Resolution: ', num2str(diff_range), ' meters']);

% Generate chirp signal for N_chirps
phase_chirp = 2 * pi * (f_start * t + (B / (2 * T_chirp)) * mod(t, T_chirp).^2);
s_t = cos(phase_chirp); % Transmitted signal

% Define target distances
d1 = 1; % Target 1 at 2.2 meters
d2 = 3; % Target 2 at 2.30 meters

disp(['Distance Between Targets: ', num2str(abs(d2 - d1)), ' meters']);

% Compute round-trip delays
tau1 = 2 * d1 / c;
tau2 = 2 * d2 / c;

% Compute sample shifts
N_shift1 = round(tau1 * Fs);
N_shift2 = round(tau2 * Fs);

% Generate received signals by shifting the transmitted chirp
s_received1 = circshift(s_t, N_shift1);
s_received2 = circshift(s_t, N_shift2);

% Superimpose the received signals
s_received = s_received1 + s_received2;
beat_signal = s_t .* s_received;

% Apply low-pass filter
beat_freq_max = range2beat(5, slope);  % Beat frequency at 5 meters
beat_signal = lowpass(beat_signal, 2 * pi * beat_freq_max, Fs);

% Plot signals
figure;
subplot(2, 1, 1)
plot(t * 1e6, s_t, 'b', 'LineWidth', 2); hold on;
plot(t * 1e6, s_received, 'r', 'LineWidth', 2);
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Transmitted\: and\: Received\: Signals}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 N_chirps * T_chirp * 1e6]) % Adjust for multiple chirps

subplot(2, 1, 2)
plot(t * 1e6, beat_signal, 'k', 'LineWidth', 2);
title("$$\bf{Mixed\: Signal}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 N_chirps * T_chirp * 1e6]) 

% Compute FFT of the beat signal
N_fft = 2^nextpow2(NumTotalPoints); % Use total points for FFT
f_axis = linspace(0, Fs / 2, N_fft / 2); % Frequency axis

% Perform FFT and normalize
S_beat_FFT = abs(fft(beat_signal, N_fft)); % Compute FFT magnitude
S_beat_FFT = S_beat_FFT(1:N_fft / 2); % Keep positive frequencies
S_beat_FFT = S_beat_FFT ./ max(S_beat_FFT); % Normalize

% Convert beat frequency to range
range_axis = beat2range(f_axis', slope);

% Plot Range Spectrum
figure;
plot(range_axis, S_beat_FFT, 'b', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Normalized\: Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Range\: Spectrum}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 5]);


% Define minimum and maximum range of interest
Rmin = 1; % Example: 0.5 meters
Rmax = 2; % Example: 4.0 meters

% Convert range limits to beat frequency using given slope
f_min = range2beat(Rmin, slope);
f_max = range2beat(Rmax, slope);

% Find indices corresponding to these frequencies
idx_min = find(f_axis >= f_min, 1, 'first'); % Index where range >= Rmin
idx_max = find(f_axis <= f_max, 1, 'last');  % Index where range <= Rmax

% Compute number of points in the selected range
num_points_in_range = idx_max - idx_min + 1;

% Output results
fprintf('Number of points between %.2f m and %.2f m: %d\n', Rmin, Rmax, num_points_in_range);

