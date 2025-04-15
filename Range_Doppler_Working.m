%% FMCW Radar Simulation: Range-Doppler Processing
clear; clc;

addpath("utils_functions_kelley");
Kelley.prettygraphs;

%This is a script to demonstrate IQ signal processing, however, it does
%not utilize IQ in it...

%% Radar Parameters
fs = 3.1e9;                % Sampling frequency (Hz)
T = 1.5e-3;                  % Chirp duration (s)
f_start = 0.5e9;            % Start frequency (Hz)
B = 1e9;                   % Bandwidth (Hz)
c = 3e8;                   % Speed of light (m/s)
sweep_slope = B / T;       % Chirp slope (Hz/s)
SNR_dB = 30;  % Signal-to-Noise Ratio in dB; set to Inf to disable noise

USE_IQ = false;  % Set to true to enable IQ processing (signed Doppler), false for standard real-valued


%% Target Parameters
R0 = 10;                   % Initial range (m)
v = -10;                     % Target velocity (m/s)

%% Simulation Parameters
num_chirps = 20;          % Number of chirps (slow-time)
t = 0:1/fs:T-1/fs;         % Fast-time vector
N = length(t);             % Samples per chirp

%% Generate Transmit Chirp 
tx = cos(2*pi * (f_start * t + 0.5 * sweep_slope * t.^2));

%% Create Mixed Signal Matrix (Mix Matrix) using circshift
mix_matrix = zeros(num_chirps, N);

for n = 1:num_chirps
    R_i = R0 + v * (n-1) * T;               % Updated range
    tau_i = 2 * R_i / c;
    delay_samples = round(tau_i * fs);
    rx_i = circshift(tx, delay_samples);

    % Add noise if SNR is finite
    if isfinite(SNR_dB)
        signal_power = rms(rx_i)^2;
        noise_power = signal_power / (10^(SNR_dB / 10));
        noise = sqrt(noise_power) * randn(size(rx_i));
        rx_i = rx_i + noise;
    end
    mix_i = tx .* rx_i;

    mix_matrix(n, :) = mix_i;
end

subplot(2,1,1)
plot(mix_matrix(1,:))
xlim([0 400])

subplot(2,1,2)
plot(mix_matrix(5,:))
xlim([0 400])

%% Optional IQ Conversion (Unambiguous Velocity via Hilbert Transform)
% Converts the real-valued signal to analytic signal if USE_IQ = true

if USE_IQ
    fprintf("Using IQ processing via Hilbert transform...\n");
    signal_matrix = hilbert(mix_matrix.').';  % Row-wise Hilbert transform
else
    fprintf("Using standard real-valued processing (Doppler ambiguity expected)...\n");
    signal_matrix = mix_matrix;
end


%% Optional: Apply windowing
win_fast = repmat(hamming(N).', num_chirps, 1);
win_slow = 1;%repmat(hamming(num_chirps), 1, N);
mix_windowed = signal_matrix .* win_fast .* win_slow;


%% Parameters for 2D FFT Padding
k_doppler = 3;  % Zero-padding factor in slow-time (Doppler) dimension
k_range   = 1;  % Zero-padding factor in fast-time (Range) dimension

padded_doppler = k_doppler * num_chirps;
padded_range   = k_range * N;

%% 2D FFT (Doppler + Range)
% Apply zero-padding in both dimensions
mix_fft2 = fftshift(fft2(mix_windowed, padded_doppler, padded_range), 1);  % fft2 and shift along slow-time

% Keep only positive range bins
mix_fft2 = mix_fft2(:, 2:floor(padded_range / 2));


%% Axes after padding
range_axis = ((0:floor(padded_range/2)) * fs / padded_range) * (c / (2 * sweep_slope));


lambda = c / (f_start + B/2);
v_max = lambda / (4 * T);  % Max unambiguous velocity
velocity_axis = linspace(-v_max, v_max, padded_doppler);


%% Plot Range-Doppler Map with Ground Truth Lines
cutoff_dB = 10;  % Set dynamic range threshold

% Compute magnitude in dB
magnitude_dB = 20 * log10(abs(mix_fft2));
magnitude_dB = max(magnitude_dB, max(magnitude_dB(:)) - cutoff_dB);  % Clip below max - cutoff

% Plot Range-Doppler Map
clf;
imagesc(range_axis, velocity_axis, magnitude_dB);
xlabel('Range (m)');
ylabel('Velocity (m/s)');
title(['Range-Doppler Map (Clipped at ', num2str(cutoff_dB), ' dB below peak)'], FontSize=25);
colormap jet;
colorbar;
axis xy;
xlim([0 2*R0]);


% --- Overlay true range and velocity ---
hold on;

% Compute expected beat frequency and Doppler shift → expected bin positions
expected_range = R0;  % initial range
expected_velocity = v;

% Draw vertical line at true range
xline(expected_range, 'w--', 'LineWidth', 1.5, 'Label', 'True Range', 'LabelVerticalAlignment', 'middle');

% Draw horizontal line at true velocity
yline(expected_velocity, 'w--', 'LineWidth', 1.5, 'Label', 'True Velocity', 'LabelHorizontalAlignment', 'center');

hold off;


