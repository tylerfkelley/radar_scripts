
addpath("utils_functions_kelley\");

Kelley.prettygraphs;
Kelley.tidy_calc;

clc; clear;

%% Parameters
c = 3e8;                    % Speed of light (m/s)
fc = 77e9;                  % Carrier frequency (Hz)
lambda = c / fc;
B = 200e6;                  % Bandwidth
T_chirp = 40e-6;            % Chirp duration
slope = B / T_chirp;
Fs = 2 * B;                 % Sampling rate
N = round(Fs * T_chirp);    % Samples per chirp
t = (0:N-1)/Fs;             % Fast time

M = 10;                     % Number of chirps (slow time)
R = 100;                    % Target range (m)
v = -10;                     % Target radial velocity (m/s)

% Compute beat frequency from range
f_beat = 2 * slope * R / c;

% Compute Doppler phase drift per chirp
delta_phi = 2 * pi * (2 * v * T_chirp * fc / c);

% Blackman window
window = blackman(N).';

% Create CPI (beat matrix)
cpi = zeros(M, N);

for m = 1:M
    noise = 0.1 * randn(1, N);
    phase_shift = (m-1) * delta_phi;
    cpi(m, :) = (sin(2*pi*f_beat*t + phase_shift) + noise) .* window;
end

%% Plot one range spectrum
rand_chirp = cpi(randi(M), :);
range_fft = fft(rand_chirp);
range_mag_dB = 20*log10(abs(range_fft) / max(abs(range_fft)));

range_axis = (0:N-1) * (c / (2 * slope * N)) * Fs;

figure;
plot(range_axis, range_mag_dB);
xlabel('Range (m)');
ylabel('Magnitude (dB)');
title('Range Spectrum (Simplified Model)');
xlim([0 200]);

%% Range-Doppler Map
range_fft = fft(cpi, [], 2);                 % Range FFT (no shift)
RD = fftshift(fft(range_fft, [], 1), 1);     % Doppler FFT (with shift)
RD_dB = 20*log10(abs(RD) / max(abs(RD(:))));

% Doppler axis
PRF = 1 / T_chirp;
fd_axis = (-M/2 : M/2 - 1) * (PRF / M);
v_axis = fd_axis * lambda / 2;

% Plot
figure;
imagesc(range_axis, v_axis, RD_dB);
xlabel('Range (m)');
ylabel('Velocity (m/s)');
title('Range-Doppler Map (Simplified Beat Model)');
colorbar;
caxis([-40 0]);
axis xy;
xlim([0 200]);

%% Peak Estimate
[max_val, idx] = max(RD_dB(:));
[row, col] = ind2sub(size(RD_dB), idx);
fprintf('Peak Detected at:\n  Range = %.2f m\n  Velocity = %.2f m/s\n', ...
        range_axis(col), v_axis(row));
