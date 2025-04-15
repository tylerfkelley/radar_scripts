% This is an FMCW simulation and also to practice OOP
clear; close all;

mat_file = 'coppertest_4G_1us_100ns.mat';

if exist(mat_file, 'file')
    load(mat_file, 'data_input'); % Load from .mat if it exists
else
    data_input = csvread("coppertest_4G_1us_100ns.csv"); % Load from CSV if .mat doesn't exist
    save(mat_file, 'data_input'); % Save to .mat for future runs
end
%%

FONTSIZE = 20;
gca_linewidth = 2.5;
set(0, 'defaultfigurecolor', [1 1 1]); % White background for plots
set(0, 'defaultAxesFontSize', FONTSIZE);
set(0, 'defaultlinelinewidth', 2);
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');

c = 3e8; % Speed of light

NumPoints = size(data_input,1);
NumScan   = size(data_input,2);

% Define start and stop frequencies of the waveform
f_start = 0.4e9; % 0.4 GHz
f_stop = 4e9; % 4 GHz
T_chirp = 1e-6; % Chirp duration (1 microsecond)
Fs = 256e9; % Sample frequency
t = (0:NumPoints - 1) / Fs; % Extended time vector
ts = 1 / Fs;

% Compute bandwidth and chirp slope
B = f_stop - f_start;
slope = B / T_chirp;
fft_resolution = 1 / T_chirp;

%% We now need to take the FFT of every row...

window_function = blackmanharris(NumPoints);
window_function = kaiser(NumPoints,5);
data_input_windowed = [];

% Use the predefined values
N_fft = 2^nextpow2(NumPoints); % Use next power of 2 for FFT efficiency
f_axis = linspace(0, Fs / 2, N_fft / 2); % Frequency axis
range_axis = beat2range(f_axis', slope); % Range axis

% Preallocate matrix for FFT results
S_beat_FFT = zeros(N_fft / 2, NumScan); % Use NumScan for clarity

% Compute FFT for each column (A-scan)
for col = 1:NumScan
    data_input_windowed(:, col) = data_input(:, col) .* window_function;
    fft_result = abs(fft(real(data_input_windowed(:, col)), N_fft)); % Compute FFT magnitude
    fft_result(1:300) = fft_result(1:300) .* 0.01;
     fft_result(1000:end) = fft_result(1000:end) .* 0.01;
     S_beat_FFT(:, col) = fft_result(1:N_fft / 2) ./ max(fft_result); % Keep positive frequencies and normalize
end

%Now i want to interpolate the FFT 


interp_factor = 10;
N_interp = size(S_beat_FFT, 1) * interp_factor; % New number of points after interpolation

% Apply interpolation column-wise
S_beat_FFT_interp = interpft(S_beat_FFT, N_interp, 1); % Interpolates along rows (frequency axis)

% Define the new interpolated frequency axis
f_axis_interp = linspace(0, Fs / 2, N_interp);

% Compute the interpolated range axis
range_axis_interp = beat2range(f_axis_interp', slope);



ascan_plotted = 60;

% Plot one of 
figure;
subplot(4, 1, 1)
hold on;
plot(t * 1e6, data_input(:,ascan_plotted), 'b', 'LineWidth', 2);
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{Time\:Domain\:Received\:Signal}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([0 0.1]) % Adjust for multiple chirps

subplot(4, 1, 2)
stem(range_axis, S_beat_FFT(:,ascan_plotted), 'b', 'LineWidth', 2);
xlabel("$$\bf{Distance\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{FFT\:of\:Received\:Signal\:(STEM)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([6 8])

subplot(4, 1, 3)
plot(range_axis, S_beat_FFT(:,ascan_plotted), 'b', 'LineWidth', 2);
xlabel("$$\bf{Distance\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{FFT\:of\:Received\:Signal\:(OLD)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([6 8])

subplot(4, 1, 4)
plot(range_axis_interp, S_beat_FFT_interp(:,ascan_plotted), 'b', 'LineWidth', 2);
xlabel("$$\bf{Distance\: (m)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
ylabel("$$\bf{Amplitude}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
title("$$\bf{FFT\:of\:Received\:Signal\:(INTERPOLATED)}$$", 'FontWeight', 'bold', 'FontSize', FONTSIZE, 'Interpreter', 'latex');
grid on;
xlim([6 8])

%% Now we start plotting some scan stuff;

x_dim = linspace(0,1.5,NumScan);

% Find indices for range between 0 and 2 meters
range_idx = find(range_axis_interp >= 6 & range_axis_interp <= 7.5);

% Extract the portion of S_beat_FFT_interp corresponding to 0-2m
S_beat_FFT_interp_trimmed = S_beat_FFT_interp(range_idx, :);
range_axis_interp_trimmed = range_axis_interp(range_idx);

S_beat_FFT_interp_trimmed = S_beat_FFT_interp_trimmed - mean(S_beat_FFT_interp_trimmed,2);

figure();
subplot(1,2,1)
plot_D(S_beat_FFT_interp_trimmed, x_dim, range_axis_interp_trimmed,"Raw Data")
xlabel("$$\bf{Scan\: x \:(m)}$$" , 'FontWeight','bold' ,FontSize=20, FontName="times");
ylabel("$$\bf{Range\: (m)}$$", 'FontWeight','bold',FontSize=20, FontName="times");
set(gca,'FontWeight','bold', FontName="times")








%% Plotting Functions Change these

function normalized_vector = normalize_to_range(input_vector, new_min, new_max)
    % Normalize a vector to a specified range [new_min, new_max]
    % 
    % Parameters:
    %   input_vector - The vector to be normalized
    %   new_min - The minimum value of the new range
    %   new_max - The maximum value of the new range
    %
    % Returns:
    %   normalized_vector - The input vector normalized to the specified range
    
    % Find the current range of the input vector
    old_min = min(input_vector);
    old_max = max(input_vector);
    
    % Scale and shift the input vector
    normalized_vector = ((input_vector - old_min) / (old_max - old_min)) * ...
                         (new_max - new_min) + new_min;
end


function x = normalize(x)
    x = x./max(abs(x(:)));
end

function [ax] = plot_I(I,x,z,DR,title_str)
colormap(jet);

imagesc(x,z.*1e2,10.*log10(abs(normalize(I))));
ax = gca();
xlabel("x (m)");
ylabel("z (cm)");
title(title_str);
caxis([-DR 0]);


end

function [ax] = plot_D(D,x,t,title_str)
    % create colorbar limits
    clim_abs_max = 250;
    threshold    = 0.99;
    D_sort = sort(D(:),'ascend');
    clims  = [D_sort(round(numel(D_sort)*(1-threshold))), D_sort(round(numel(D_sort)*threshold))];
    clims  = [max(-clim_abs_max, clims(1)),min(clim_abs_max, clims(2))];

    colormap(gray(2^12));
    imagesc(x, t, D);
    ax = gca();
    %colorbar;
    caxis(clims);
    title(title_str, FontSize=20);
    xlabel("x (m)", FontSize=20);
    ylabel("f (Ghz)", FontSize=20);
end