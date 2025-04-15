%This is an fmcw sim and also to practice OOP
clear; close all;


FONTSIZE = 20;
gca_linewidth = 2.5;
set(0,'defaultfigurecolor',[1 1 1]); %white background for plots...
set(0,'defaultAxesFontSize',FONTSIZE); 
set(0,'defaultlinelinewidth',2);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');


c = 3e8;

Fs = 256e9; %Sample freq.
NumPoints = 512e3; %Number of points

t = (0:NumPoints-1) / Fs; % Time vector

% Define start and stop frequencies
f_start = 0.4e9; % 0.4 GHz
f_stop = 4e9;  % 4 GHz
T_chirp = 1e-6; % Chirp duration of 1 microsecond


% Compute bandwidth
B = f_stop - f_start;
slope = B / T_chirp;

% Convert range limits to beat frequencies
beat_freq_min = range2beat(0.1, slope);  % Beat frequency at 0 meters
beat_freq_max = range2beat(5, slope);  % Beat frequency at 2 meters

% Generate chirp signal
phase_chirp = 2 * pi * (f_start * t + (B / (2 * T_chirp)) * t.^2); % Instantaneous phase


s_t = cos(phase_chirp); % FMCW chirp signal


% Define target distances
d1 = 2.2; % Target 1 at 0.2 meters
d2 = 2.2; % Target 2 at 0.5 meters

% Compute round-trip delays
tau1 = 2 * d1 / c;
tau2 = 2 * d2 / c;

% Compute sample shifts
N_shift1 = round(tau1 * Fs);
N_shift2 = round(tau2 * Fs);

% Generate received signals by shifting the transmitted chirp
s_received1 = circshift(s_t, N_shift1); % Signal from target 1
s_received2 = circshift(s_t, N_shift2); % Signal from target 2

% Superimpose the received signals
s_received = s_received1 + s_received2;
beat_signal = s_t .* s_received;

beat_signal = lowpass(beat_signal,2*pi*beat_freq_max,Fs);

figure;
subplot(2,1,1)
plot(t * 1e6, s_t, 'b', 'LineWidth', 2); hold on;
plot(t * 1e6, s_received, 'r', 'LineWidth', 2);
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
%legend('Transmitted Chirp', 'Received Signal (Two Targets)', 'Location', 'best');
title("$$\bf{Theoretical\: Transmitted\: and\: Received\: Signals}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
xlim([0 0.01])
set(gca, 'LineWidth', gca_linewidth, 'FontWeight', 'bold', 'FontName', 'times');
subplot(2,1,2)
plot(t * 1e6, beat_signal, 'k', 'LineWidth', 2);
title("$$\bf{Theoretical\: Mixed\: Signal}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
xlim([0 0.01])


figure;
set(gca, 'LineWidth', gca_linewidth, 'FontWeight', 'bold', 'FontName', 'times');
plot(t * 1e6, beat_signal, 'k', 'LineWidth', 2);
title("$$\bf{Theoretical\: Mixed\: Signal}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;


% Compute FFT of the beat signal
N_fft = 2^nextpow2(NumPoints); % Zero-padding for better resolution
f_axis = linspace(0, Fs/2, N_fft/2); % Frequency axis (positive frequencies)

% Perform FFT and normalize
S_beat_FFT = abs(fft(beat_signal, N_fft)); % Compute FFT magnitude
S_beat_FFT = S_beat_FFT(1:N_fft/2); % Keep only positive frequencies
S_beat_FFT = S_beat_FFT ./ max(S_beat_FFT);

% Convert beat frequency to range
range_axis = beat2range(f_axis', slope); % Convert frequencies to range

% Plot Range Spectrum
figure;
plot(range_axis, S_beat_FFT, 'b', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Normalize Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
title("$$\bf{Theoretical\: Range\: Spectrum}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
xlim([0 5])


%% This is to read some real data to do an experiment

input = csvread("cabletest_channel2_4G_2_25_2025.csv");
t = input(:,1);
input = input(:,2);

output = csvread("cabletest_function1_4G_2_25_2025.csv");
output = output(:,2);


ts = mean(diff(t));
fs = 1/ts;


% Define start and stop frequencies
f_start = 0.4e9; % 0.4 GHz
f_stop = 4e9;  % 4 GHz
T_chirp = 1e-6; % Chirp duration of 1 microsecond


% Compute bandwidth
B = f_stop - f_start;
slope = B / T_chirp;

% Convert range limits to beat frequencies
beat_freq_min = range2beat(0.1, slope);  % Beat frequency at 0 meters
beat_freq_max = range2beat(3, slope);  % Beat frequency at 2 meters

mixed = input .* output;

 %mixed = highpass(mixed,2*pi*beat_freq_min,fs);
 %mixed = lowpass(mixed,2*pi*beat_freq_max,fs);

figure;
subplot(2,1,1)
plot(t * 1e6, input, 'b', 'LineWidth', 0.5); hold on;
plot(t * 1e6, output, 'r', 'LineWidth', 0.5);
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
legend('Transmitted Chirp', 'Received Signal (Two Targets)', 'Location', 'best');
title("$$\bf{Experimental\: Transmitted\: and\: Received\: Signals}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
set(gca, 'LineWidth', gca_linewidth, 'FontWeight', 'bold', 'FontName', 'times');
subplot(2,1,2)
plot(t * 1e6, mixed, 'k', 'LineWidth', 2);
title("$$\bf{Experimental\: Mixed\: Signal}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;

% Compute FFT of the beat signal
N_fft = 2^nextpow2(length(t)); % Zero-padding for better resolution
f_axis = linspace(0, fs/2, N_fft/2); % Frequency axis (positive frequencies)

% Perform FFT and normalize
Mixed = abs((fft(mixed, N_fft))); % Compute FFT magnitude
Mixed = Mixed(1:N_fft/2); % Keep only positive frequencies
%Mixed = 10*log10(Mixed);

Mixed = Mixed ./ (max(Mixed));

% Convert beat frequency to range
range_axis = beat2range(f_axis', slope); % Convert frequencies to range

% Plot Beat Frequency Spectrum
figure;
plot(range_axis, Mixed, 'k', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Normalize Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
title("$$\bf{Experimental\: Range\: Spectrum}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
xlim([0 5])

figure();
plot(f_axis./1e6, Mixed, 'k', 'LineWidth', 2);
xlabel("$$\bf{Frequency\: (MHz)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Normalize Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
title("$$\bf{Experimental\: Range\: Spectrum}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
xlim([0 100])


%% This is to read some real data to do an experiment

input = csvread("boarsight_channel2_4G_2_25_2025.csv");
t = input(:,1);
input = input(:,2);
% input  = diff(input);
% input(end+1) = 0;

output = csvread("boarsight_function1_4G_2_25_2025.csv");
output = output(:,2);



ts = mean(diff(t));
fs = 1/ts;

% Define start and stop frequencies
f_start = 0.4e9; % 0.4 GHz
f_stop = 4e9;  % 4 GHz
T_chirp = 1e-6; % Chirp duration of 1 microsecond


% Compute bandwidth
B = f_stop - f_start;
slope = B / T_chirp;

% Convert range limits to beat frequencies
beat_freq_min = range2beat(0.1, slope);  % Beat frequency at 0 meters
beat_freq_max = range2beat(3, slope);  % Beat frequency at 2 meters

mixed = input .* output;

% mixed = highpass(mixed,2*pi*beat_freq_min,fs);
% mixed = lowpass(mixed,2*pi*beat_freq_max,fs);

figure;
subplot(2,1,1)
plot(t * 1e6, input, 'b', 'LineWidth', 0.5); hold on;
plot(t * 1e6, output, 'r', 'LineWidth', 0.5);
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
legend('Transmitted Chirp', 'Received Signal (Two Targets)', 'Location', 'best');
title("$$\bf{Boresight\: Transmitted\: and\: Received\: Signals}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
set(gca, 'LineWidth', gca_linewidth, 'FontWeight', 'bold', 'FontName', 'times');
subplot(2,1,2)
plot(t * 1e6, mixed, 'k', 'LineWidth', 2);
title("$$\bf{Boresight\: Mixed\: Signal}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;

% Compute FFT of the beat signal
N_fft = 2^nextpow2(length(t)); % Zero-padding for better resolution
f_axis = linspace(0, fs/2, N_fft/2); % Frequency axis (positive frequencies)

% Perform FFT and normalize
Mixed = abs((fft(mixed, N_fft))); % Compute FFT magnitude
Mixed = Mixed(1:N_fft/2); % Keep only positive frequencies
%Mixed = 10*log10(Mixed);


Mixed = Mixed ./ (max(Mixed));
%Mixed = 10*log10(Mixed);

% Convert beat frequency to range
range_axis = beat2range(f_axis', slope); % Convert frequencies to range

% Plot Beat Frequency Spectrum
figure;
plot(range_axis, Mixed, 'k', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Normalize Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
title("$$\bf{Boresight\: Range\: Spectrum}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
xlim([4 7])

%% This is when S-Parameters are applied:

%% This is to read some real data to do an experiment

input = csvread("cabletest_channel2_4G_2_25_2025.csv");
t = input(:,1);
input = input(:,2);

output = csvread("cabletest_function1_4G_2_25_2025.csv");
output = output(:,2);


ts = mean(diff(t));
fs = 1/ts;

output = apply_sparam_to_waveform('multilayer_1meter.s2p',output,ts);

% Define start and stop frequencies
f_start = 0.4e9; % 0.4 GHz
f_stop = 4e9;  % 4 GHz
T_chirp = 1e-6; % Chirp duration of 1 microsecond


% Compute bandwidth
B = f_stop - f_start;
slope = B / T_chirp;

% Convert range limits to beat frequencies
beat_freq_min = range2beat(0.1, slope);  % Beat frequency at 0 meters
beat_freq_max = range2beat(3, slope);  % Beat frequency at 2 meters

mixed = input .* output;

 %mixed = highpass(mixed,2*pi*beat_freq_min,fs);
 %mixed = lowpass(mixed,2*pi*beat_freq_max,fs);


% Compute FFT of the beat signal
N_fft = 2^nextpow2(length(t)); % Zero-padding for better resolution
f_axis = linspace(0, fs/2, N_fft/2); % Frequency axis (positive frequencies)

% Perform FFT and normalize
Mixed = abs((fft(mixed, N_fft))); % Compute FFT magnitude
Mixed = Mixed(1:N_fft/2); % Keep only positive frequencies
%Mixed = 10*log10(Mixed);

Mixed = Mixed ./ (max(Mixed));

% Convert beat frequency to range
range_axis = beat2range(f_axis', slope); % Convert frequencies to range;

figure();
plot(f_axis./1e6, Mixed, 'k', 'LineWidth', 2);
xlabel("$$\bf{Frequency\: (MHz)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Normalize Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
title("$$\bf{S-Parameters\: From\: Cable\: Test}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
xlim([0 100])

%% Do this again but try out some different frequencies
c = 3e8
Fs = 256e9; %Sample freq.
NumPoints = 256e3; %Number of points

t = (0:NumPoints-1) / Fs; % Time vector

% Define start and stop frequencies
f_start = 4e9; % 0.4 GHz
f_stop = 7e9;  % 4 GHz
T_chirp = 1e-6; % Chirp duration of 1 microsecond


% Compute bandwidth
B = f_stop - f_start;
slope = B / T_chirp;

% Convert range limits to beat frequencies
beat_freq_min = range2beat(0.1, slope);  % Beat frequency at 0 meters
beat_freq_max = range2beat(5, slope);  % Beat frequency at 2 meters

% Generate chirp signal
phase_chirp = 2 * pi * (f_start * t + (B / (2 * T_chirp)) * t.^2); % Instantaneous phase


s_t = cos(phase_chirp); % FMCW chirp signal


% Define target distances
d1 = 2.2; % Target 1 at 0.2 meters
d2 = 2.2; % Target 2 at 0.5 meters

% Compute round-trip delays
tau1 = 2 * d1 / c;
tau2 = 2 * d2 / c;

% Compute sample shifts
N_shift1 = round(tau1 * Fs);
N_shift2 = round(tau2 * Fs);

% Generate received signals by shifting the transmitted chirp
s_received1 = circshift(s_t, N_shift1); % Signal from target 1
s_received2 = circshift(s_t, N_shift2); % Signal from target 2

% Superimpose the received signals
s_received = s_received1 + s_received2;

s_received = apply_sparam_to_waveform('multilayer_1meter.s2p',s_received',1/Fs);


beat_signal = s_t .* s_received';

beat_signal = lowpass(beat_signal,2*pi*beat_freq_max,Fs);


% Compute FFT of the beat signal
N_fft = 2^nextpow2(NumPoints); % Zero-padding for better resolution
f_axis = linspace(0, Fs/2, N_fft/2); % Frequency axis (positive frequencies)

% Perform FFT and normalize
S_beat_FFT = abs(fft(beat_signal, N_fft)); % Compute FFT magnitude
S_beat_FFT = S_beat_FFT(1:N_fft/2); % Keep only positive frequencies
S_beat_FFT = S_beat_FFT ./ max(S_beat_FFT);

% Convert beat frequency to range
range_axis = beat2range(f_axis', slope); % Convert frequencies to range

% Plot Range Spectrum
figure;
plot(range_axis, S_beat_FFT, 'b', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Normalize Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
title("$$\bf{Theoretical\: Range\: Spectrum}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
xlim([0 5])


%% This is to read some real data to do an experiment

input = csvread("boarsight_channel2_2.5G_2_25_2025.csv");
t = input(:,1);
input = input(:,2);
% input  = diff(input);
% input(end+1) = 0;

output = csvread("boarsight_function1_2.5G_2_25_2025.csv");
output = output(:,2);



ts = mean(diff(t));
fs = 1/ts;

% Define start and stop frequencies
f_start = 0.4e9; % 0.4 GHz
f_stop =  4e9;  % 4 GHz
T_chirp = 1e-6; % Chirp duration of 1 microsecond


% Compute bandwidth
B = f_stop - f_start;
slope = B / T_chirp;

% Convert range limits to beat frequencies
beat_freq_min = range2beat(4, slope);  % Beat frequency at 0 meters
beat_freq_max = range2beat(7, slope);  % Beat frequency at 2 meters

mixed = input .* output;



figure;
subplot(2,1,1)
plot(t * 1e6, input, 'b', 'LineWidth', 0.5); hold on;
plot(t * 1e6, output, 'r', 'LineWidth', 0.5);
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
legend('Transmitted Chirp', 'Received Signal (Two Targets)', 'Location', 'best');
title("$$\bf{Boresight\: Transmitted\: and\: Received\: Signals}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
set(gca, 'LineWidth', gca_linewidth, 'FontWeight', 'bold', 'FontName', 'times');
subplot(2,1,2)
plot(t * 1e6, mixed, 'k', 'LineWidth', 2);
title("$$\bf{Boresight\: Mixed\: Signal}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;

% Compute FFT of the beat signal
N_fft = 2^nextpow2(length(t)); % Zero-padding for better resolution
f_axis = linspace(0, fs/2, N_fft/2); % Frequency axis (positive frequencies)

% Perform FFT and normalize
Mixed = abs((fft(mixed, N_fft))); % Compute FFT magnitude
Mixed = Mixed(1:N_fft/2); % Keep only positive frequencies
%Mixed = 10*log10(Mixed);


% Normalize the input signal
Mixed = abs(Mixed) ./ max(abs(Mixed)); % Normalize between 0 and 1

% CFAR Parameters
N_train = 3;  % Number of training cells (reduce if signal disappears)
N_guard = 1;   % Number of guard cells (reduce if needed)
P_fa = 1e-1;   % Desired false alarm probability
alpha = (N_train * (P_fa^(-1/N_train)) - 1); % Adaptive threshold scaling

% Initialize CFAR output
Mixed_cfar = zeros(size(Mixed));

% CFAR Sliding Window Process
for i = (N_train + N_guard + 1):(length(Mixed) - (N_train + N_guard))
    % Training region (excluding guard cells)
    training_cells = [Mixed(i - N_train - N_guard : i - N_guard - 1); ...
                      Mixed(i + N_guard + 1 : i + N_train + N_guard)];
    
    % Estimate noise level using median (more robust than mean)
    noise_power = median(training_cells);
    
    % Compute CFAR threshold
    threshold = alpha * noise_power;
    
    % Apply CFAR detection
    if Mixed(i) > threshold
        Mixed_cfar(i) = Mixed(i); % Keep detected target
    end
end

% Plot CFAR results
figure;
plot(Mixed, 'b', 'LineWidth', 1.2); hold on;
plot(Mixed_cfar, 'r', 'LineWidth', 2);
xlabel('Frequency Bin');
ylabel('Normalized Amplitude');
title('CFAR Detection on Frequency Spectrum');
legend('Original Spectrum', 'Detected Targets');
grid on;





%Mixed = 10*log10(Mixed);

% Convert beat frequency to range
range_axis = beat2range(f_axis', slope); % Convert frequencies to range

% Plot Beat Frequency Spectrum
figure;
plot(range_axis, Mixed_cfar, 'k', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Normalize Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
title("$$\bf{Boresight\: Range\: Spectrum}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
xlim([4 7])



%% Lets plot the s-paramters and the phase.

sparam_file = "multilayer_1meter.s2p";

S = sparameters(sparam_file); % Read S-parameters from file

% Extract frequency, S21 magnitude, and phase
freq = S.Frequencies; % Frequency in Hz
S21 = squeeze(S.Parameters(2,1,:)); % Extract S21 (2nd port, 1st port)

% Convert magnitude to dB
S21_dB = 20*log10(abs(S21)); % Magnitude in dB

% Extract phase in degrees and unwrap it
S21_phase = angle(S21) * (180/pi); % Convert from radians to degrees
S21_phase_unwrapped = unwrap(S21_phase * pi/180) * (180/pi); % Unwrap and convert back to degrees

% Compute group velocity
delta_f = diff(freq); % Frequency spacing (Hz)
delta_phi = diff(S21_phase_unwrapped) * (pi/180); % Convert degrees to radians
group_velocity = (2 * pi * delta_f) ./ delta_phi; % Compute group velocity

% Compute midpoint frequency for plotting
freq_mid = (freq(1:end-1) + freq(2:end)) / 2;

% Plot S21 magnitude
figure;
subplot(2,1,1);
plot(freq/1e9, S21_dB, 'b', 'LineWidth', 2);
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('S21 Magnitude vs Frequency');
grid on;


% Plot group velocity
subplot(2,1,2);
plot(freq_mid/1e9, group_velocity, 'r', 'LineWidth', 2);
xlabel('Frequency (GHz)');
ylabel('Group Velocity (m/s)');
title('Group Velocity vs Frequency');
grid on;

%% Lets plot the s-paramters from the system

sparam_file = 'cleaned_system_fmcw.s2p';
S = sparameters(sparam_file); % Read S-parameters from file

% Extract frequency, S21 magnitude, and phase
freq = S.Frequencies; % Frequency in Hz
S21 = squeeze(S.Parameters(2,1,:)); % Extract S21 (2nd port, 1st port)

% Convert magnitude to dB
S21_dB = 20*log10(abs(S21)); % Magnitude in dB

% Extract and unwrap phase (keep in radians)
S21_phase_unwrapped = unwrap(angle(S21)); % Unwrap phase (radians)

% Compute group velocity
delta_f = diff(freq); % Frequency spacing (Hz)
delta_phi = diff(S21_phase_unwrapped); % Phase difference (radians)
group_velocity = - (2 * pi * delta_f) ./ delta_phi; % Compute group velocity
group_velocity = log10(group_velocity);
% Compute midpoint frequency for plotting
freq_mid = (freq(1:end-1) + freq(2:end)) / 2;

% Plot S21 magnitude
figure;
subplot(2,1,1);
plot(freq/1e9, S21_dB, 'b', 'LineWidth', 2);
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('S21 Magnitude vs Frequency');
grid on;

% Normalize if needed (optional)
group_velocity = group_velocity ./ max(abs(group_velocity)); 

% Plot group velocity
subplot(2,1,2);
plot(freq_mid/1e9, group_velocity, 'r', 'LineWidth', 2);
xlabel('Frequency (GHz)');
ylabel('Normalized Group Velocity');
title('Group Velocity vs Frequency');
grid on;

%%

%% Do this again but try out some different frequencies
c = 3e8
Fs = 256e9; %Sample freq.
NumPoints = 256e3; %Number of points

t = (0:NumPoints-1) / Fs; % Time vector

% Define start and stop frequencies
f_start = 0.4e9; % 0.4 GHz
f_stop = 5e9;  % 4 GHz
T_chirp = 1e-6; % Chirp duration of 1 microsecond


% Compute bandwidth
B = f_stop - f_start;
slope = B / T_chirp;

% Convert range limits to beat frequencies
beat_freq_min = range2beat(0.1, slope);  % Beat frequency at 0 meters
beat_freq_max = range2beat(5, slope);  % Beat frequency at 2 meters

% Generate chirp signal
phase_chirp = 2 * pi * (f_start * t + (B / (2 * T_chirp)) * t.^2); % Instantaneous phase


s_t = cos(phase_chirp); % FMCW chirp signal




% Define target distances
d1 = 2.2; % Target 1 at 0.2 meters
d2 = 2.2; % Target 2 at 0.5 meters

% Compute round-trip delays
tau1 = 2 * d1 / c;
tau2 = 2 * d2 / c;

% Compute sample shifts
N_shift1 = round(tau1 * Fs);
N_shift2 = round(tau2 * Fs);

% Generate received signals by shifting the transmitted chirp
s_received1 = circshift(s_t, N_shift1); % Signal from target 1
s_received2 = circshift(s_t, N_shift2); % Signal from target 2

% Superimpose the received signals
s_received = s_received1 + s_received2;

s_received = apply_sparam_to_waveform('cleaned_system_fmcw.s2p',s_received',1/Fs);


beat_signal = s_t .* s_received';

%beat_signal = lowpass(beat_signal,2*pi*beat_freq_max,Fs);


% Compute FFT of the beat signal
N_fft = 2^nextpow2(NumPoints); % Zero-padding for better resolution
f_axis = linspace(0, Fs/2, N_fft/2); % Frequency axis (positive frequencies)

% Perform FFT and normalize
S_beat_FFT = abs(fft(beat_signal, N_fft)); % Compute FFT magnitude
S_beat_FFT = S_beat_FFT(1:N_fft/2); % Keep only positive frequencies
S_beat_FFT = S_beat_FFT ./ max(S_beat_FFT);

% Convert beat frequency to range
range_axis = beat2range(f_axis', slope); % Convert frequencies to range

% Plot Range Spectrum
figure;
plot(range_axis, S_beat_FFT, 'b', 'LineWidth', 2);
xlabel("$$\bf{Range\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
ylabel("$$\bf{Normalize Magnitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
title("$$\bf{Theoretical\: Range\: Spectrum}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'FontName', 'times', 'Interpreter', 'latex');
grid on;
xlim([0 5])

%%


function output_waveform = apply_sparam_to_waveform(sparam_file, your_waveform, Ts)
% APPLY_SPARAM_TO_WAVEFORM Applies S-parameter transfer function to a waveform
%
% Inputs:
% sparam_file - String, path to the .s2p or .s1p file containing S-parameters
% your_waveform - Vector, time-domain signal to be processed
% Ts - Scalar, sampling period of your waveform
%
% Output:
% output_waveform - Vector, processed waveform after applying S-parameter response
% Load S-parameters
S = sparameters(sparam_file);
S_matrix = S.Parameters; % Extract S-parameter matrix
freq = S.Frequencies; % Frequency points
% Extract S21 (for a two-port system, assumes transmission response)
H = squeeze(S_matrix(2,1,:));
% Define waveform frequency axis
N = length(your_waveform);
fs_waveform = 1/Ts; % Sampling rate of the waveform
freq_waveform = (0:N-1)*(fs_waveform/N); % Proper frequency vector for the waveform
% Only use frequencies up to Nyquist (fs/2)
nyquist_freq = fs_waveform/2;
valid_freq_indices = freq_waveform <= nyquist_freq;
freq_waveform_valid = freq_waveform(valid_freq_indices);
% Find valid S-parameter frequencies (within our signal bandwidth)
valid_sparam_indices = freq <= nyquist_freq;
freq_valid = freq(valid_sparam_indices);
H_valid = H(valid_sparam_indices);
% Handle case where S-parameters don't start from DC (0 Hz)
if freq_valid(1) > 0
    freq_valid = [0; freq_valid];
    H_valid = [H_valid(1); H_valid]; % Assume DC response is same as lowest frequency
end
% Interpolate S21 to match waveform frequency points
H_interp = interp1(freq_valid, H_valid, freq_waveform_valid, 'linear', 'extrap');
% Create complete frequency response (including negative frequencies for ifft)
H_full = zeros(N, 1);
H_full(valid_freq_indices) = H_interp;
% Handle negative frequencies (complex conjugate of positive frequencies)
% Skip DC (index 1) and Nyquist (if N is even)
if mod(N, 2) == 0 % even number of points
    neg_indices = N:-1:N/2+2;
    pos_indices = 2:N/2;
else % odd number of points
    neg_indices = N:-1:(N+3)/2;
    pos_indices = 2:(N+1)/2;
end
H_full(neg_indices) = conj(H_full(pos_indices));
% Compute FFT of the waveform
Y_waveform = fft(your_waveform);
% Apply the frequency response
Y_out = Y_waveform .* H_full;
% Convert back to time domain
output_waveform = real(ifft(Y_out));
end















% function output_waveform = apply_sparam_to_waveform(sparam_file, your_waveform, Ts)
%     % APPLY_SPARAM_TO_WAVEFORM Applies S-parameter transfer function to a waveform
%     % 
%     % Inputs:
%     %   sparam_file   - String, path to the .s2p or .s1p file containing S-parameters
%     %   your_waveform - Vector, time-domain signal to be processed
%     %   Ts            - Scalar, sampling period of your waveform
%     % 
%     % Output:
%     %   output_waveform - Vector, processed waveform after applying S-parameter response
%     
%     % Load S-parameters
%     S = sparameters(sparam_file);
%     S_matrix = S.Parameters; % Extract S-parameter matrix
%     freq = S.Frequencies; % Frequency points
%     
%     % Extract S21 (for a two-port system, assumes transmission response)
%     H = squeeze(S_matrix(2,1,:));
%     
%     % Define waveform frequency axis
%     fs_waveform = 1/Ts; % Sampling rate of the waveform
%     N = length(your_waveform);
%     freq_interp = linspace(min(freq), max(freq), N);
%     
%     % Interpolate S21 to match waveform frequency points
%     H_interp = interp1(freq, H, freq_interp, 'linear', 'extrap');
%     
%     % Compute FFT of the waveform
%     Y_waveform = fft(your_waveform);
%     
%     % Apply the frequency response
%     Y_out = Y_waveform .* H_interp(:); % Ensure correct dimensions
%     
%     % Convert back to time domain
%     output_waveform = ifft(Y_out, 'symmetric');
%     
%  
% end


