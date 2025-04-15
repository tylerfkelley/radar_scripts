% FmcwSim.m

% Follows conventions in "Fourier Transform and Its Applications" by R.N. Bracewell
% Note: f_x denotes a function of the real space variable x.
% Note: X_w denotes a function of the angular frequency variable w.
% 10/25/2024

% If the figure already exists clear it, otherwise create it
try
	if exist('fig_s', 'var') 
		clf(fig_s);			
	else
		fig_s = figure();
	end
catch
	clear all;
	close all;

	fig_s = figure();
end

% Create the OneDimensionalTransform object and setup the independent variables
% Real space extent, number of points (try even/odd number of points to include zero)
UsPerChirp = 10;		% This is the real space extent, the real space values will range from +/- 1/2 this value
NumPoints = 1001;	% Use an odd number here
odt = OneDimensionalTransforms(UsPerChirp, NumPoints);	

StartFrequencyMHz = 1;
ChirpRateMhzPerUs = 2/10;

% Example of chirp usage. Inputs are starting number of cycles per waveform, number of cycles to chirp, shift
% [f_x, NumShiftSamples] = odt.Chirp(CyclesPerChirp, DeltaCyclesPerChirp, 0); 
[f_x, NumShiftSamples] = odt.Chirp(StartFrequencyMHz, ChirpRateMhzPerUs, 0); 

[X_w, f_x_ht] = odt.HilbertTransform(f_x, 0);
f_x_ht = (f_x_ht);

figure(fig_s);
NumPlots = 5;
PlotType = 'line'; % 'stem'
% Plot the chirp
subplot(NumPlots,1,1)
odt.PlotReImX(f_x, PlotType, "real(f(x))", "imag(f(x))",  '$f(x)$',  '$Time \ (\mu S)$');
hold on
% Plot the hilbert transform of the chirp
subplot(NumPlots,1,2)
odt.PlotReImX(f_x_ht, PlotType, "real(f(x))", "imag(f(x))",  '$f(x)$',  '$Time \ (\mu S)$');
plot(odt.x, abs(f_x_ht));
legend(["real(HT(f(x)))", "imag(HT(f(x)))", "Envelope(f(x)"]);
grid(gca,'minor')
grid on;


% Create the delayed response
DelayUs = 5;
DelayShift = floor(DelayUs*(NumPoints/UsPerChirp));
f_x_rx = circshift(f_x, DelayShift);
f_x_rx(1:DelayShift) =0;

% Mix the tx and rx
f_x_if = f_x_rx .* f_x;

% Plot the transmitted and recieved chirp as well as the mixed IF
subplot(NumPlots,1,3)
plot(odt.x, f_x);
hold on
plot(odt.x, f_x_rx);
plot(odt.x, f_x_if);
xlabel('$Time \ (\mu S)$')
legend(["TX", "RX", "IF"])

% Filter the mixed response
PadFactor = 2;
[X_w_if, Window] = odt.FourierTransform(f_x_if, PadFactor);


% Filter the response
MaxFilterBin = 60;
FftLength = length(X_w_if);
X_w_if(MaxFilterBin:(FftLength-MaxFilterBin)) = 0;

% Transform back to real space
[f_x_if] = odt.InverseFourierTransform(X_w_if, Window, PadFactor);


subplot(NumPlots,1,4)
plot(odt.x, real(f_x_if))
xlabel('$Time \ (\mu S)$')
legend(["Filtered IF"])



% Keep the first part of the IF specturm
NumIfPoints = floor(NumPoints/10);
IfSpectrum = abs(X_w_if(1:NumIfPoints));
IfFrequencies = odt.FrequencyInterp(1:NumIfPoints);
IfTimeDelaysUs = IfFrequencies/ChirpRateMhzPerUs;

MetersPerUs = 299.792;
RangeInMeters = (1/2)*IfTimeDelaysUs*MetersPerUs;

[PeakY, PeakX] = max(IfSpectrum);
IfFrequency = IfFrequencies(PeakX);
TimeOffsetUs = IfFrequency/ChirpRateMhzPerUs;



subplot(NumPlots,1,5)
% plot(IfFrequencies, IfSpectrum)
% xlabel("IF Frequency (Mhz)")
plot(IfTimeDelaysUs, IfSpectrum)
xlabel("IF TimeOffset (Us)")
% plot(RangeInMeters, IfSpectrum)
% xlabel("Range (m)")

% odt.PlotMagPhaseW(X_w_if, PlotType, "$\phi(X(\omega))$", "$|X(\omega)|$", '$\phi(X(\omega)) \ and \ (X(\omega))|$', 'Frequency');
% xlim([min(IfFrequencies) max(IfFrequencies)]);







