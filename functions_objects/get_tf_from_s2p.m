function H = get_tf_from_s2p(filename, f_desired)
% Samuel Wagner, UC Davis ECE MML, 2021
% this function will take an .s2p's filename as an input and
% return its transfer function over a desired frequency range
%
% this likely requires RF toolbox!
%
% inputs
% filename - filename 
% f_desired - the frequency vector with which to conform H
% outputs
% H - voltage transfer function

% read the s2p and get its frequencies
S = sparameters(filename);
f = S.Frequencies;

% throw some warnings if frequency ranges don't add up nicely.
if(max(f) < max(f_desired))
    warning("Requested frequency of " + num2str(max(f_desired)./1e9) + ...
        " GHz is higher than max file frequency of " + num2str(max(f)./1e9) + ...
        " GHz. Zeroing out non-prescribed frequencies using an FIR-20 filter.");
end
if(min(f)  > min(abs(f_desired)))
    warning("Requested frequency of " + num2str(min(abs(f_desired))./1e6) + ...
        " MHz is lower than file min frequency of " + num2str(min(f)./1e6) + ...
        " MHz. Zeroing out non-prescribed frequencies using an FIR-20 filter.");
end

% calculate transfer function using s2tf
tf 	= s2tf(S);

% tf = squeeze(S.Parameters(2,1,:));
% reshape vectors, ensure stability of interpolation
f 	= f(:);
tf 	= tf(:);
tf 	= tf(f>0);
f   = f(f>0);

% reshape f and tf s.t. they are hermetian
f_H  = [-flip(f); f];
tf_H = [flip(flip_phase(tf)); tf];

%interpolate tf over f_desired
H = interp1(f_H, tf_H, f_desired, 'linear','extrap');

% unset the non-trusted frequencies
if(max(f) < max(f_desired))
    Filter_b = fir1(50, max(f)./max(f_desired).*0.90);
    TF = freqz(Filter_b,1, f_desired, max(f_desired)*2);
    H = H.*TF;
end
