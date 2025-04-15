
%Tyler Kelley
%This is an oject to run simple FMCW simulations and test theory.

classdef FMCWsimulation < handle
    
    properties

    NumPoints; %Number of sample points in a chirp.
    
    Fs;        %Sample frequency

    x;         %Spatial vector

    beat;     %Beat frequency results

    chirp;    %FMCW transmitted signal

    mixed;    %FMCW mixed signal

    target;   %Return signal offset by a certain distance

    t;        %time

    end



%Create an initialization for the FMCW parameters
    methods
        function obj = FMCWsimulation(num_points, fs)
            % Constructor: Initializes parameters
            if nargin > 0
                obj.NumPoints = num_points;
                obj.Fs = fs;
            else
                obj.NumPoints = 1024;   % Default number of points
                obj.Fs = 1e9;           % Default sampling frequency (1 GHz)
            end
            
            % Create the time vector
            obj.x = linspace(0, obj.NumPoints / obj.Fs, obj.NumPoints);
        end


        function obj = generateChirp()





    end
end