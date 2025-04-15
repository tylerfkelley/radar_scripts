% Brendan Illingworth
% Follows conventions in "Fourier Transform and Its Applications" by R.N. Bracewell
% Note: f_x denotes a function of the real space variable x.
% Note: X_w denotes a function of the angular frequency variable w.
% 12/24/2022


classdef OneDimensionalTransforms < handle

properties

	SpatialExtent;			% Real space variable ranges from -SpatialExtent/2 to +SpatialExtent/2
	NumPoints;				% Number of points in the independent variables
	
	x;						% Real space independent variable

	DftW;					% Fourier space independent variable in DFT terms e.g. [0, 2pi]
	Omega;					% Fourier space independent variable in angular frequency
	Frequency;				% Fourier space independent variable in frequency units
	FrequencyInterp;		% Fourier space independent variable in frequency units after zero padding

	RealStep;				% Distance between samples in real space			
	DftStepW;				% Distance between samples in DFT space	in radians
	DftStepF;				% Distance between samples in DFT space	in frequency units

	SampleRate;				% SampleRate in real space
	DeltaF;					% Frequency domain resolution in Hz

	MaxAngularFrequency		% Max angular frequency
	DeltaOmega;				% Frequency domain resolution in radians

	CenterIndex;

end
   
methods

	% Create the real space and Fourier space independent variables, for some functions it is convenient to exclude zero
	function obj = OneDimensionalTransforms(SpatialExtent, NumPoints)

		% Set default linewidth to draw thicker lines
		set(0, 'DefaultLineLineWidth', 2);
		
		% Use latex interpreter by default
		set(groot,'defaultAxesTickLabelInterpreter','latex');  
		set(groot,'defaulttextinterpreter','latex');
		set(groot,'defaultLegendInterpreter','latex');

		% Setup the indepenedent variables
		obj.SpatialExtent = SpatialExtent;
		obj.NumPoints = NumPoints;
		obj.RealStep = SpatialExtent/(NumPoints-1);
		obj.DftStepW = 2*pi/NumPoints;
		obj.DftStepF = 1/SpatialExtent;

		% Define the time and frequency axis in terms of actual values
		obj.SampleRate = 1/(obj.RealStep);
		obj.DeltaF = obj.SampleRate/(obj.NumPoints);			 
		obj.DeltaOmega = 2*pi*obj.DeltaF;
		obj.MaxAngularFrequency = obj.DeltaOmega*(NumPoints/2);		% Half the points are the negative frequencies
 

		% Ensure that we don't include zero for reason above 
		obj.x = linspace(-obj.SpatialExtent/2, obj.SpatialExtent/2, obj.NumPoints);

		% Frequency space independent variable in frequency units 
		obj.Frequency = obj.DftStepF*(0:1:(obj.NumPoints-1));

		% Ensure that in frequency space omega ranges from 0 to 2*pi minus the frequency step size! 
		obj.DftW = linspace(0, (2*pi)-obj.DftStepW, obj.NumPoints);
		
		% Create angular frequency vector positive
		obj.Omega = 0:obj.DeltaOmega:obj.MaxAngularFrequency;

		% Create angular frequency vector with positive/negative frequencies
		obj.Omega = [obj.Omega -flip(obj.Omega(2:end))];					
		obj.Omega = fftshift(obj.Omega);						% 02/27 this convention requires rotating FFT array

		obj.CenterIndex = floor(obj.NumPoints/2) + 1;
	end

%% Basic Functions
	% Return one of the basic functions stretched based on width
	function [y, NumSamples] = GetBasicFunction(obj, func, width, shift)

		switch func

			case 'rect'
				y = obj.Rect(width);
			case 'triangle'
				y = obj.Triangle(width);
			case 'gaussian'
				y = obj.Gaussian(width);
			case 'GaussianTau'
				y = obj.GaussianTau(width);
			case 'cos'
				y = obj.Cos(width);
			case 'CosOmega'
				y = obj.CosOmega(width);
			case 'sin'
				y = obj.Sin(width);
			case 'sync'
				y = obj.Sinc(width);
			case 'SechTau'
				y = obj.SechTau(width);	
			case 'comb'
				% The fundamental period of the comb function is extactly one delta function per spatial extent, with 'harmonics' equally
				% spaced. The comb function requires the scale be such we evaluate on integer mutliples of the spatial step size.
				a = obj.RealStep*10;
				y = obj.Comb((width + obj.RealStep)*a);
			otherwise

		end

		% Now shift by the correct number of samples
		NumSamples = floor(shift/obj.RealStep);
  		y = circshift(y, NumSamples);

	end

	% The following functions are naturally aperiodic and the first input represents the width: 'rect', 'triangle', 'gaussian', 'sync'
	% The following functions are naturally periodic and the first input respresents the number of cycles: 'cos', 'sin'
	% The following functions are naturally periodic and the first input represents the number of sample between periods: 'comb'

	% Rectangle function of unit height and base (non-zero between 0 and 1)
	% Use shift and stretch to spatially shift and stretch spatial width
    function y = Rect(obj, width)

		% Length of a row is the number of collumns 
		y = zeros(1, length(obj.x));

		% The unstretched rect function has this many non-zero samples between -1/2 and +1/2
		NumSamples = 1/obj.RealStep;
		NumSamples = floor(NumSamples * width);		% scale

		% Create rect between 0 and 1
		y(obj.CenterIndex:obj.CenterIndex+NumSamples) = 1;

	end

	% Triangle function of unit height and base (non-zero between -1/2 and +1/2)
	% Use shift and stretch to spatially shift and stretch spatial width
	function y = Triangle(obj, width)

		% Length of a row is the number of collumns 
		y = zeros(1, obj.NumPoints);

		% The unstretched triangle function has this many non-zero samples between -1/2 and +1/2
		NumSamples = floor(width/obj.RealStep);

		% Create the triangle centered around zero (the center of the spatial axis)
		indices = obj.CenterIndex-floor(NumSamples/2):1:obj.CenterIndex+floor(NumSamples/2);
		y(indices) = 2*(1/2-abs(obj.x(indices)/width));

	end

	% Zero shift puts the function in the center of the array
	function y = Gaussian(obj, width)
		y = exp(-pi*(width*obj.x).^2);
	end

	% Zero shift puts the function in the center of the array
	function y = GaussianTau(obj, tau)

		t = obj.x;
		y = exp(-(t.^2)./(2*(tau.^2)));

	end

	% Note that the starting phase will depending on the spatial axis setup. 
	function y = Cos(obj, NumCyclesInWaveform)
		% Scale the angular frequency so that we get the desired number of cycles per waveform
		omega = 2*pi*NumCyclesInWaveform/(obj.SpatialExtent);
		y = cos(omega*obj.x);
	end

	% Note that the starting phase will depending on the spatial axis setup. 
	function y = Sin(obj, NumCyclesInWaveform)
		% Scale the angular frequency so that we get the desired number of cycles per waveform
		omega = 2*pi*NumCyclesInWaveform/(obj.SpatialExtent);
		y = sin(omega*obj.x);
	end

	% Zero shift puts the function in the center of the array
	function y = Sinc(obj, width)
		y = sinc(width*obj.x);
	end

	% Zero shift puts the function in the center of the array
	function y = SechTau(obj, tau)
		y = sech(obj.x/tau);
	end

	% Comb(a(x))
	function y = Comb(obj, width)

		% sum over n of delta((a(x))-n)
		indicies = find(mod(obj.x, width) == 0);

		% Length of a row is the number of collumns 
		y = zeros(1, length(obj.x));
		y(indicies) = 1;

	end

%% More Complex Functions

	function [y, NumSamples]= Chirp(obj, StartFrequencyMHz, ChirpRateMhzPerUs, shift)

		UsPerChirp = obj.SpatialExtent;
		NumCyclesStart = StartFrequencyMHz * UsPerChirp;

		StopFrequencyMhz = StartFrequencyMHz + ChirpRateMhzPerUs * UsPerChirp;
		NumCyclesChirp = (StopFrequencyMhz-StartFrequencyMHz) * UsPerChirp;

		% Scale the angular frequency so that we get the desired number of cycles per waveform
		omega = 2*pi*NumCyclesStart/(obj.SpatialExtent);
		
		% Create a chirp to add to the start frequency
		phase = linspace(0, sqrt(NumCyclesChirp*pi), obj.NumPoints);
		y = cos(omega*obj.x + phase.^2);

		% Now shift by the correct number of samples
		NumSamples = floor(shift/obj.RealStep);
  		y = circshift(y, NumSamples);

	end

	% Zero pad, shift, window and transform
	function [X_w, Window]= FourierTransform(obj, f_x, ZeroPadFactor)

		% Length of zero padding and amount of shift to center signal in the window
		PadLength = ZeroPadFactor*obj.NumPoints;
		shift = floor(PadLength/2);

		% Zero padding increases the frequency resolution 
		obj.FrequencyInterp = (obj.DftStepF/(1+ZeroPadFactor))*(0:1:(obj.NumPoints+PadLength-1));


		% Create the window function
		% Window = [];
		Window = blackmanharris(obj.NumPoints + PadLength)';

		% 1) Zero pad 
		pad = zeros(1, PadLength);
		x = [f_x pad];

		% 2) Shift so the function in the center of the window
		x = circshift(x, shift);
		
		% 3) Apply the window
		if Window
			x = x.*Window;
		else
			x = x;
		end

		% 4) Take the Fourier Transform
		X_w = fft(x);
		
	end

	% Inverse transform, invert window, shift
	function [f_x_ft] = InverseFourierTransform(obj, X_w, Window, ZeroPadFactor)

		% Length of zero padding and amount of shift to center signal in the window
		PadLength = ZeroPadFactor*obj.NumPoints;
		shift = floor(PadLength/2);

		% 4) Inverse transfrom back to real space
		f_x_ft = ifft(X_w);

		% 3) Invert the window
		if Window
			f_x_ft = f_x_ft ./ Window;
		end

		% 2) Shift the real space function back
		f_x_ft = circshift(f_x_ft, -shift);

		% 1) Remove the zero padding
		f_x_ft = f_x_ft(1:obj.NumPoints);

		
	end

	function [X_w, f_x_ht] = HilbertTransform(obj, f_x, ZeroPadFactor)

		[X_w, Window] = FourierTransform(obj, f_x, ZeroPadFactor);
		
		% Find the new center index
		ci = (1+ZeroPadFactor)*obj.CenterIndex;

		% Zero all the negative frequencies
		X_w(1:ci) = 0;

		% Double the magnitude of the positive frequencies to preserve power
		X_w(ci:end) = 2*X_w(ci:end);

		f_x_ht = InverseFourierTransform(obj, X_w, Window, ZeroPadFactor);
		
	end

%% Ploting Functions

	% Plot real/imag part of spatial function
	function PlotReImX(obj, f_x, type, legend1, legend2, Title, xTitle)

		if strcmp(type, 'stem')
			stem(obj.x, real(f_x), 'filled')
			hold on;
			stem(obj.x, imag(f_x),'filled')
		else
			plot(obj.x, real(f_x))
			hold on;
			plot(obj.x, imag(f_x))
		end

 		title(Title);
		legend([legend1, legend2]);
		ylabel('Amplitude');
		xlabel(xTitle);
		xlim([obj.x(1), obj.x(end)])

	end

	% Plot real/image part of angular frequency function
	function PlotReImW(obj, X_w, type, legend1, legend2, Title)

		if strcmp(type, 'stem')
			stem(obj.DftW, real(X_w), 'filled')
			hold on;
			stem(obj.DftW, imag(X_w),'filled')
		else
			plot(obj.DftW, real(X_w))
			hold on;
			plot(obj.DftW, imag(X_w))
		end

 		title(Title);
		legend([legend1, legend2]);
		ylabel('Amplitude');
		xlabel('$\omega$');

		% Replace the default numberical axis text with latex symbols, note that last element is 2*pi-stepsize
		xlim([obj.DftW(1), obj.DftW(end)])
		xticks([0 pi/2 pi 3*pi/2 obj.DftW(end)])
		xticklabels({'0', '$\frac{\pi}{2}$', '${\pi}$','$\frac{3\pi}{2}$','$2{\pi}$'})
	end

	% Plot magnitude/phase of angular frequency function
	function PlotMagPhaseW(obj, X_w, type, legend1, legend2, Title, xType)

		if strcmp(xType, 'Omega')
			xvals = obj.DftW;
		else		
			if length(X_w) == length(obj.Frequency)
				xvals = obj.Frequency;
			else
				xvals = obj.FrequencyInterp;
			end
		end

		if strcmp(type, 'stem')
			yyaxis left
			stem(xvals, unwrap(angle(X_w)), 'filled')
			ylabel('$Phase \ in \ Radians$')
			yyaxis right		
			stem(xvals, abs(X_w), 'filled')
		else
			yyaxis left
			plot(xvals, unwrap(angle(X_w)))
			ylabel('$Phase \ in \ Radians$')
			yyaxis right		
			plot(xvals, abs(X_w))
		end

 		title(Title);
		legend([legend1, legend2]);
		ylabel('Amplitude');


		if strcmp(xType, 'Omega')
			xlabel('$\omega$');
			% Replace the default numberical axis text with latex symbols, note that last element is 2*pi-stepsize
			xlim([xvals(1), xvals(end)])
			xticks([0 pi/2 pi 3*pi/2 xvals(end)])
			xticklabels({'0', '$\frac{\pi}{2}$', '${\pi}$','$\frac{3\pi}{2}$','$2{\pi}$'})
		else
			xlabel(xType);
		end

	end

	%
	function delete(obj)

	end

end

end

			









