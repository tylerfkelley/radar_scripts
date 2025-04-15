%%
% Script NAME: Kelley_Interview_Question.m
%
% AUTHOR: Kelley, Tyler
% DATE:   2 / 20 / 2025
% PREPARED FOR: Anduril Interview Question on radar design.

% Problem set:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a concept for an FMCW radar that meets the requirements stated below.
% Write an analysis script  to predict the range at which the seeker will detect the spec target.
% Start by showing a plot of Signal-to-Noise ratio vs. range and
% compute the target acquisition range by assuming Swerling 1 target fluctuation statistics and 
% a 90% single look probability of detection.
% Include a description of a candidate waveform that will be used for target acquisition.
% State any assumptions needed to complete the analysis.
% These can include system parameters such as Tx power, Receiver noise figure, 
% and Transmit-Receive antenna isolation.
%
% Spec target RCS: 0.1 meter^2
% Min waveform Range: 10 meters
% Max waveform Range: 200 meters
% Range resolution: 1 meter
% Max target speed: 300 m/s
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% This is for the plots to make everything look nice.

clear; close all;


FONTSIZE = 20;
gca_linewidth = 2.5;
set(0,'defaultfigurecolor',[1 1 1]); %white background for plots...
set(0,'defaultAxesFontSize',FONTSIZE); 
set(0,'defaultlinelinewidth',2);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

%% Variable Initilization

%These are first the variables given in the problem set.
   sigma            = 0.1; %meter^2 (RCS)
   min_range        = 10;  %meter
   max_range        = 200; %meter  
   range_resolution = 1;   %meter
   max_speed        = 300; %meters per second
   Pd               = 0.9;  %Probability of detection
   k                = 1.38e-23; % Boltzmann s Constant (J/K)
   T                = 290; % Standard Temperature (Kelvin)
   c                = 3e8; % Speed of light (m/s)

 %Select propagation environment condition
 propagation_condition = 'clear'; % Options: 'clear', 'rain', 'fog'

%We now have to make some assumptions about our RF seeker to continue on
%with the calculation.

tidy_calc;
disp(['---------Start of Script Processing---------']);
tidy_calc;
disp(['---------Assumed Parameters---------']);

%We will first assume the center frequency of the seeker. As referenced:
%https://apps.dtic.mil/sti/tr/pdf/ADA317256.pdf
%The KU band (12 - 18 GHz) can be implemented:

fc = 15e9; % Center frequency (Hz)
lambda = c / fc; % Wavelength (m)

%We then need to assume the transmit power. 

Pt = 10; %Transmitted power in watts;

%Based off of the following refernce:
%https://apps.dtic.mil/sti/tr/pdf/ADA111647.pdf
%We can get a realistic assumption for antenna gain based on the seeker
%size.

aperature_efficiency = 0.90;
seeker_diameter      = 20e-2; %meters
effective_antenna_aperature = aperature_efficiency*pi*(seeker_diameter/2)^2;

G = (4*pi*effective_antenna_aperature) ./ lambda^2;
G_dB = lin2db(G);
disp(['Antenna gain: ', num2str(G_dB), ' dBi']);

%Based off of the following refernce:
%https://apps.dtic.mil/sti/tr/pdf/ADA111647.pdf
% as well as sample questions from "Principles of Modern Radar"
%We can make some further assumptions for a typical receiver noise figure
%as well as system losses...

%These are assumed values for noise figures in the receiving chain.
% Noise Figures of components (in dB)

NF_LIM   = 1.0;   % Limiter (assumed 2Vpp limiter)
NF_LNA   = 1.5;   % LNA
NF_BPF   = 0.5;   % BPF
NF_MIXER = 3.0;   % Mixer
NF_ADC   = 10.0;  % ADC (Usually the highest NF)

% Gains of each stage (in dB)
G_LIM   = -1;   % Limiter Loss (insertion loss)
G_LNA   = 20;   % LNA Gain
G_BPF   = -1.5; % BPF Loss
G_MIXER = -6;   % Mixer Loss
G_ADC   =  0;   % ADC

% Convert Noise Figures to Noise Factors (Linear Scale)
F_LIM   = 10^(NF_LIM/10);
F_LNA   = 10^(NF_LNA/10);
F_BPF   = 10^(NF_BPF/10);
F_MIXER = 10^(NF_MIXER/10);
F_ADC   = 10^(NF_ADC/10);

% Convert Gains to Linear Scale
G_LIM   = 10^(G_LIM/10);
G_LNA   = 10^(G_LNA/10);
G_BPF   = 10^(G_BPF/10);
G_MIXER = 10^(G_MIXER/10);
G_ADC   = 10^(G_ADC/10);

% Compute Total Gain
G_total = G_LIM * G_LNA * G_MIXER * G_BPF * G_ADC;

% Compute Total Noise Factor using Friis Formula
Nf = F_LNA + ...
     (F_MIXER - 1) / G_LNA + ...
     (F_BPF - 1) / (G_LNA * G_MIXER) + ...
     (F_LIM - 1) / (G_LNA * G_MIXER * G_BPF) + ...
     (F_ADC - 1) / (G_LNA * G_MIXER * G_BPF * G_LIM);


Nf_dB = db2lin(Nf); %Noise Figure in dB

%Give a generic number for system losses.
L_dB  = 2;  %System Losses in dB
L  = db2lin(L_dB);

disp(['Noise Figure: ', num2str(Nf_dB), ' dB']);
disp(['System Losses: ', num2str(L_dB), ' dB']);


%We can assume some additional losses from Tx / Rx isolation:

Tx_Rx_Isolation_dB = 60; % Isolation in dB
Tx_Rx_Isolation = db2lin(Tx_Rx_Isolation_dB);
disp(['Tx / Rx Isolation: ', num2str(Tx_Rx_Isolation_dB), ' dB']);

%We also need to assume the false alarm probability. Based off the
%following:
%https://123.physics.ucdavis.edu/week_5_files/filters/matched_filter.pdf
%We choose the typical value of:

P_fa = 1e-6;
disp(['False Alarm Probability: ', num2str(P_fa), '']);

%Values of 10^-4 and 10^-8 also seem to be common.

%A lower Pfa will give fewer false alarms but also requires a higher SNR. 
%Essentially there is a tradeoff.

%We also need to assume this device will not just operate in free space.

%We first need to define the range vector.
range_vector = linspace(min_range,max_range,1000);
%1000 points is chosen arbitrarily...


% Compute atmospheric attenuation based on condition
switch propagation_condition
    case 'clear'
        L_atm_dB = 2 * gaspl(range_vector, fc, 15, 1013*1e2, 7.5); % T = 15°C
        disp('Assuming propagation losses in air on a clear day.');

    case 'rain'
        L_atm_dB = 2 * rainpl(range_vector, fc, 95); % Heavy rain: 95 mm/hr
        disp('Assuming heavy rain conditions (95 mm/hr rainfall).');

    case 'fog'
        % Fog attenuation using MATLAB's fogpl() function
        L_atm_dB = 2 * fogpl(range_vector, fc, 15, 0.8); % T = 15°C, Liquid water density = 0.8 g/m³
        disp('Assuming dense fog conditions (0.8 g/m³ liquid water).');

    otherwise
        error('Invalid propagation condition. Choose "clear", "rain", or "fog".');
end




L_atm = db2lin(L_atm_dB);



%% These are variables for the waveform we will have to design based on the
%problem set.

tidy_calc;
disp(['---------Calculated Waveform Parameters---------']);

%This is decided from the basic range resolution equation.
B = c / (2 * range_resolution); % Bandwidth (Hz)
disp(['Required Bandwidth: ', num2str(B/1e6), ' MHz']);

%This first step is to decide essentially the minimum sweep time needed to
%avoid range ambiguity.

T_sweep_min = (2 * max_range) / c;
sweep_factor = 5.5;
T_sweep_min = T_sweep_min * sweep_factor;
disp(['Minimum Sweep Time (Range limited) : ', num2str(T_sweep_min*1e6), ' µs']);

%Is referenced in: https://www.mathworks.com/help/radar/ug/automotive-adaptive-cruise-control-using-fmcw-technology.html
%That the minium sweep time should have atleast an additional factor of 5.5
% to ensure the beat frequency is within receivers sampling range,
% preventing aliasing.


%There is also a doppler contraint with the sweep time we need to consider
%If we made the sweep time really long, we couldn't detect the targets
%at 300 m/s.

T_sweep_max = lambda / (4*max_speed);
disp(['Maximum Sweep Time (Velocity Limited) : ', num2str(T_sweep_max*1e6), ' µs']);

%I now choose the sweep time to be about 90% of the max so we don't miss
%high velocity targets but also improve matched filter gain.

T_sweep = T_sweep_max .* 0.9;


disp(['Chosen Sweep Time : ', num2str(T_sweep*1e6), ' µs']);

disp(['Chosen Modulation Scheme : Triangular']);



%% Now I will start by plotting the signal-to-noise ratio against range.

%Assumptions for signal processing and range gating are also made here.

%First define the thermal noise power coupled with reciever noise figure:
Pn = k * T * Nf * B;

% Calculate coupling power (Antenna crosstalk)
Pc = Pt / Tx_Rx_Isolation; % Coupling power in watts

tidy_calc;
disp(['---------Crosstalk Assumptions---------']);


%This is admittedly where I spent the most time with this problem and where
%I could rationalize a few different methods of how to model the cross talk
%or the antenna isolation in this problem. In thie first case, we could
%model any of the cross talk as noise since it would essentially be this
%constant signal the reciever would be dealing with.
%In the second case, we could also model it as noise added to the system,
%however, I can assume we would also have the ability to place an analog
%filter after the mixer in the system to get rid of any lower beat
%frequencies. Essentially, since the cross talk would just be seen as a
%really strong close target, and our minumum range is 10 meters, I would
%just design a highpass filter to get rid of any below 10 meter beat
%frequencies which would attenuate that "noise". Finally, and I think the
%way that I would ultimately decide, is to use the analog filter but not
%treat the antenna coupling as noise in the system. Instead, what I
%ultimately think that matters most is thinking about the ADC in the
%system. A limiter could be very helpful.

 crosstalk_condition = 'coupling loss'; % Options: 'noise', 'range gated', 'coupling loss'

switch crosstalk_condition
    case 'noise'
        %In this case, we simply add in the antenna cross talk to the total
        %noise of the system
        Pn = Pn + Pc;
        disp(['Assume cross talk directly correlates to noise.']);
    case 'range gated'

        %Essentially we could include an analog filter to exclude some of the lower beat 
        %frequencies that would occur below 10m, mainly the coupling.
        %Consider triangular modulation.
        f_low  = range2beat(10,(2*B)/T_sweep,c);
        disp(['Assume Filtering below: ', sprintf('%.0f', f_low*1e-3), ' KHz due to fb.']);

        filter_suppression_dB = 30; %Assume 30dB supression from the filter below the cutoff.
        filter_suppression_linear = 10^(-filter_suppression_dB/10);  % Convert dB to linear scale
        Pc_filtered = Pc * filter_suppression_linear;  % Apply filtering attenuation
        Pn = Pn + Pc_filtered;  % Add filtered leakage to noise floor
        %https://www.mdpi.com/1424-8220/22/12/4641
    case 'coupling loss'

        %In this case, we consider mainly the coupling as a loss factor,
        %but ultimately have to consider the ADC in future steps such as
        %dynamic range, sampling rate, NF, ETC... 

        Pt = Pt - Pc;    
        disp(['Assume cross talk directly correlates to loss of transmit power.']);

        otherwise
        Pn;
end

%However, we also have to think about the ADC in the system before moving
%on. We wouldnt want to saturate our receiver with a potential close
%target. This is why a 2Vpp limiter is chosen in the chain. I imagine for
%an RF seeker there is a chance for the device to get very close to a
%target. Specs for the ADC will be mentioned later.



%% Plot SNR vs. Range

SNR_linear = (Pt * G^2 * lambda^2 * sigma) ./ ((4 * pi)^3 * range_vector.^4 * Pn * L .* L_atm');
SNR_dB = lin2db(SNR_linear); % Convert SNR to dB

figure(1);
plot(range_vector, SNR_dB, 'b', 'LineWidth', 3);
grid on;
title("$$\bf{SNR\: vs.\: Range\: for\: FMCW\: Seeker\:}$$" , 'FontWeight','bold' ,FontSize=22, FontName="times");
xlabel("$$\bf{Range\: (m)}$$" , 'FontWeight','bold' ,FontSize=22, FontName="times");
ylabel("$$\bf{SNR\: (dB)}$$", 'FontWeight','bold',FontSize=22, FontName="times");
set(gca,'LineWidth',gca_linewidth);
set(gca,'FontWeight','bold', FontName="times")
xlim([min_range max_range]);



%% Swerling 1 Considerations

% The Swerling models describe different levels of target fluctuations
% due to varying radar cross-section (RCS) characteristics. This affects
% how a target reflects radar signals and thus impacts detection.
% This was my first time really deep diving into the subject so these were
% my basic findings:
% Swerling 0: Constant RCS - Non-Fluctuating Target (GPR Targets)
% Swerling 1: Rayleigh Distributed RCS - each pulse sees a different RCS,
% but remains constant during a detection period.
% Swerling 2: Fast fluctuating, RCS changes between pulses.
% Since Swerling 1 targets fluctuate in RCS, their detection follows an 
% exponential distribution. The probability of detection is given by:
%
%    Pd = integral(  Q1(sqrt(2 * SNR * x), sqrt(2 * gamma)) * exp(-x) dx )
%
% This integral has no closed-form solution, so we need to solve it
% numerically. Admittedly I spent some time on this and the values were not
% seeming to work out correctly. It wasnt until thoroughly examining
% chapter 3 of "Principles of Modern Radar" I found that these equations
% are not very strightforward to solve and lookup tables are normally used.
% I then deep more searching and found that MATLAB has a built in equation
% for the Schnidman function, which is exactly what I needed in the
% situation. Shnidman's equation is a series of equations that yield an 
% estimate of the SNR required for a specified false-alarm and detection probability.
% After further research I found I could also use Albersheims equation, but
% Schnidman's was better versed for moving targets.

%%
tidy_calc;
disp(['---------Swerling 1 PDF Considerations---------']);

% Set parameters
Pd_range = 0.1:0.01:0.99;  % Probability of Detection range
SNR_dB_plotted = [];
% Calculate required SNR for Swerling 1 model
%This equation takes in the probability of detection, probability of false
%alarm rate, mumber of integreated pulses, and the Swerling model number.

%To further investigate this, I wanted to plot the required SNR over Pd
%also to verify my results.

for i = 1:length(Pd_range)
SNR_dB_plotted(i) = shnidman(Pd_range(i),P_fa,1,1);
end

% Plot the curve
fig = figure(2);
fig.Position(3) = 900;  % Change width (increase for wider figure)
plot(Pd_range, SNR_dB_plotted, 'r', 'LineWidth', 3)
grid on;
xlim([0.1 1])
xticks(0:0.1:1)
title("$$\bf{Required\: SNR\: vs.\: Pd\: for\: Swerling\: 1\: Model}$$" , 'FontWeight','bold' ,FontSize=22, FontName="times");
xlabel("$$\bf{Probability\: of\: Detection}$$" , 'FontWeight','bold' ,FontSize=22, FontName="times");
ylabel("$$\bf{Required\: SNR\: (dB)}$$", 'FontWeight','bold',FontSize=22, FontName="times");
set(gca,'LineWidth',gca_linewidth);
set(gca,'FontWeight','bold', FontName="times")
xline(Pd, 'k--', 'LineWidth', 3);
text(0.5, 10, "$$P_{fa} = 10^{-6}$$", 'FontSize', 22, 'Color', 'k', 'FontWeight', 'bold');



% Compute required SNR for Swerling 1 using statistical model.
% For a single look we essentially only coherently integrate one pulse.
SNR_req_dB = shnidman(Pd,P_fa,1,1);

disp(['Required SNR for Swerling 1 Target: ', num2str(SNR_req_dB), ' dB']);

%%

% Determine Maximum Detection Range

% Find the last range value where computed SNR meets or exceeds required SNR
detection_range = range_vector(find(SNR_dB >= SNR_req_dB, 1, 'last'));

% Display the result
if detection_range > 0
disp(['Maximum Detection Range (Swerling 1): ', num2str(detection_range), ' meters']);
else
disp(['SNR is too low!']);
end

% Plot SNR vs. Required SNR to Show Detection Range

fig = figure(3);
fig.Position(3) = 900;  % Adjust width for better visualization
plot(range_vector, SNR_dB, 'b', 'LineWidth', 3); hold on;
yline(SNR_req_dB, 'r--', 'LineWidth', 3); % Required SNR line

grid on;
xlabel("$$\bf{Range\: (m)}$$" , 'FontWeight','bold' ,'FontSize',22, 'FontName',"times");
ylabel("$$\bf{SNR\: (dB)}$$", 'FontWeight','bold', 'FontSize',22, 'FontName',"times");
title("$$\bf{SNR\: vs.\: Range\: for\: FMCW\: Seeker}$$" , 'FontWeight','bold' ,'FontSize',22, 'FontName',"times");
set(gca,'LineWidth',gca_linewidth);
set(gca,'FontWeight','bold', FontName="times")
legend('Computed SNR', 'Required SNR (Swerling 1)', 'Location', 'NorthEast');
xlim([min_range max_range]);

% Annotate the detection range
if detection_range >= 0
text(detection_range-50, SNR_req_dB + 10, ...
    strcat("Detection Range: ", num2str(detection_range), " m"), ...
    'FontSize', 16, 'Color', 'r', 'FontWeight', 'bold');
 % Display warning message if SNR is too low
else
    text(max(range_vector)/2, max(SNR_dB)/2, ...  % Position in middle of plot
        "SNR is too low!", ...
        'FontSize', 18, 'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

%% FMCW Waveform Design Explanation
%
% The radar system employs an FMCW waveform,
% which enables simultaneous range and velocity measurement. The waveform consists
% of a triangular chirped signal with:
%
% - Bandwidth:  B = 150 MHz → Determines range resolution
% - Sweep Time: T_sweep = 12 µs → Affects maximum velocity detection
%
% Why a Triangular FMCW Waveform?
% - The frequency sweeps up and then down within each cycle.
% - This allows direct measurement of both range and Doppler shift.
% - By comparing beat frequencies from up and down sweeps, the system resolves range and velocity.
% - However, this also makes us consider the chirp rate is cut in half
% - (Higher ADC sampling)

tidy_calc;
disp(['---------Necessary Hardware Specs---------']);

%https://www.mdpi.com/1424-8220/22/12/4641

%Need to consider triangular modulation so sweep time is halved.
f_low  = range2beat(min_range,B/(T_sweep/2),c);
f_high_range  = range2beat(max_range,B/(T_sweep/2),c);

%We have to also consider the maximum doppler shift.
fd_max = speed2dop(2*max_speed,lambda);

f_high = f_high_range + fd_max;

disp(['Bandpass filter needed between: ', sprintf('%.2f', f_low*1e-6),  ...
    ' and ' sprintf('%.0f', f_high*1e-6), ' MHz.']);


%This is for ADC specs

% Limiter Output (Max Power in Watts)
P_max = 0.01; % Derived from 2Vpp limit

% Compute P_min from radar equation (200 meters away)
P_min = (Pt * G^2 * lambda^2 * sigma * G_total) ./ ((4 * pi)^3 * 200.^4 .* L .* L_atm(200));

% Compute Dynamic Range (dB)
DR = 10 * log10(P_max ./ P_min);
disp(['Minimum ADC Dynamic Range: ', num2str(DR), ' dB']);

ADC_sampling_minimum = f_high .* 2.5;
disp(['Minimum ADC Sampling Rate: ', sprintf('%.0f', ADC_sampling_minimum.*1e-6), ' MSPS (2.5x Nyquist)']);


%% Waveform visualization:

% Waveform parameters
f_start = fc - B/2; % Start frequency (Hz) = 14.925 GHz
f_end = fc + B/2;   % End frequency (Hz) = 15.075 GHz
T_half = T_sweep / 2; % Half sweep time (s)
rate_half = B / T_half; % Chirp rate (Hz/s)

% Time vector and frequency
fs = ADC_sampling_minimum; % Sampling frequency (1 GHz)
t = 0:1/fs:T_sweep; % Time vector (s)
f_tri = zeros(size(t));
f_tri(t <= T_half) = f_start + rate_half * t(t <= T_half); % Up-sweep
f_tri(t > T_half) = f_end - rate_half * (t(t > T_half) - T_half); % Down-sweep

% Convert frequency to GHz for plotting
f_ghz = f_tri / 1e9; % Frequency in GHz
t_us = t * 1e6; % Time in microseconds

% Plot
figure(4);
plot(t_us, f_ghz, 'b', 'LineWidth', 3);
grid on;
title("$$\bf{Triangular\: FMCW\: Waveform}$$", 'FontWeight', 'bold', 'FontSize', 22, 'FontName', 'times');
xlabel("$$\bf{Time\: (\mu s)}$$", 'FontWeight', 'bold', 'FontSize', 22, 'FontName', 'times');
ylabel("$$\bf{Frequency\: (GHz)}$$", 'FontWeight', 'bold', 'FontSize', 22, 'FontName', 'times');
set(gca, 'LineWidth', gca_linewidth);
set(gca, 'FontWeight', 'bold', 'FontName', 'times');
set(gca, 'FontSize', FONTSIZE);
xlim([0 T_sweep*1e6]); % Full sweep time in µs
ylim([14.925 15.075]); % Frequency range in GHz


%% Functions live down here...

function linear_value = db2lin(db_value)
    linear_value = 10.^(db_value / 10);
end

function db_value = lin2db(linear_value)
    db_value = 10 * log10(linear_value);
end

function tidy_calc = tidy_calc()
%This is just to give some space in the command window when needed.
disp(['           ']);
disp(['           ']);
end












