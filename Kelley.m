classdef Kelley
    methods (Static)
        % Normalize data
        function x = normalize(x)
            x = x ./ max(abs(x(:)));
        end

        % Calculate Power Spectral Density (PSD) of a waveform
        function X = calculate_psd(x, Fs)
            X = fftshift(fft(x, [], 1), 1);
            X = 1 / (Fs * size(x, 1)) .* abs(X).^2;
        end

        % Convert Vpp to dBm (Assuming 50 Ohms)
        function dBm = convert_dbm(Vpp)
            Vp = Vpp / 2;
            Vrms = Vp * 0.707;
            P_mw = Vrms.^2 * 10^3 / 50; 
            dBm = 10 * log10(P_mw);
        end

        % Read CSV file from ADS
        function waveform = read_ads_csvfile(title, size)
            main = readtable(title);
            main = main.Variables;
            main = strrep(main, ' ', '');
            main = strrep(main, 'psec', 'e-12');
            main = strrep(main, 'fsec', 'e-15');
            main = strrep(main, 'nsec', 'e-9');
            main = strrep(main, 'aV', 'e-18');
            main = strrep(main, 'fV', 'e-15');
            main = strrep(main, 'pV', 'e-12');
            main = strrep(main, 'mV', 'e-3');
            main = strrep(main, 'nV', 'e-9');
            main = strrep(main, 'uV', 'e-6');
            main = strrep(main, 'sec', 'e0');
            main = strrep(main, 'V', 'e0');
            waveform_init = str2double(main);

            new_time_index = linspace(1, numel(waveform_init(:,1)), size);
            resized_time = interp1(waveform_init(:,1), new_time_index);
            new_data_index = linspace(1, numel(waveform_init(:,2)), size);
            resized_data = interp1(waveform_init(:,2), new_time_index);

            waveform = [resized_time', resized_data'];
        end

        % Find Full Width at Half Maximum (FWHM)
        function pulse_width = fwhm(t, x)
            if length(t) == length(x)
                x = abs(x);
                [~, index] = max(x);
                left_index = find(x(1:index) < (x(index) / 2), 1, 'last');
                right_index = find(x(index:end) < (x(index) / 2), 1, 'first') + index - 1;
                pulse_width = (t(right_index) - t(left_index)) * 1e12;
            else
                pulse_width = 0;
                warning("Both vectors must be the same size. First input is time, second is voltage.");
            end
        end

        % Extract time vector from waveform
        function x = time_vec(y)
            if size(y, 2) > 2
                error('This function is intended for two-column lists.');
            end
            x = y(:,1);
        end

        % Extract data vector from waveform
        function x = data_vec(y)
            if size(y, 2) > 2
                error('This function is intended for two-column lists.');
            end
            x = y(:,2);
        end

        % Compute PSD and normalize for plotting
        function X = psd_data(x, Fs)
            X = fftshift(fft(x, [], 1), 1);
            X = 1 / (Fs * size(x, 1)) .* abs(X).^2;
            X = 10 .* log10(Kelley.normalize(abs(X)));
            [max_X, ~] = max(X);
            X = X + abs(max_X);
        end

        % Create TIM file
        function makeTimfile(title, t, y)
            length_tim = length(t);
            first_string = "BEGIN TIMEDATA";
            sec_string = "# T ( SEC V R xx )";
            third_string = "% time voltage";
            last_string = "END";
            init_array = [t', y'];
            
            fid = fopen(title + ".tim", 'w');
            fprintf(fid, '%s\n', first_string);
            fprintf(fid, '%s\n', sec_string);
            fprintf(fid, '%s\n', third_string);
            for i = 1:length_tim
                fprintf(fid, '%s\n', string(init_array(i, 1)) + "	" + string(init_array(i, 2)));
            end
            fprintf(fid, '%s\n', last_string);
            fclose(fid);
        end

        function plot_D(D,x,t,title_str)
        % create colorbar limits
        clim_abs_max = 250;
        threshold    = 0.99;
        D_sort = sort(D(:),'ascend');
        clims  = [D_sort(round(numel(D_sort)*(1-threshold))), D_sort(round(numel(D_sort)*threshold))];
        clims  = [max(-clim_abs_max, clims(1)),min(clim_abs_max, clims(2))];
    
        figure(); colormap(gray(2^12));
        imagesc(x, t.*1e9, D);
        colorbar;
        caxis(clims);
        title(title_str);
        xlabel("x (m)");
        ylabel("t (ns)");
        end

        function plot_I(I,x,z,DR,title_str)
            figure(); colormap(jet);
            imagesc(x,z.*1e2,10.*log10(abs(normalize(I))));
            ylabel("z (cm)");
            xlabel("x (m)");
            title(title_str);
            caxis([-DR 0]);
            colorbar;
        end

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

        function X = prettygraphs()
            FONTSIZE = 20;
            set(0,'defaultfigurecolor',[1 1 1]); %white bkgr
            set(0,'defaultAxesFontSize',FONTSIZE); 
            set(0,'defaultlinelinewidth',2);
            set(0,'defaultLegendInterpreter','latex');
            set(0,'defaultTextInterpreter','latex');
        end

        function x_padded = zpad(x, multiplier, mode)
        %ZERO_PAD_BY_MULTIPLIER Zero-pads a vector to a multiple of its length.
        %   x_padded = ZERO_PAD_BY_MULTIPLIER(x, multiplier, mode)
        %
        %   x         : Input vector (row or column)
        %   multiplier: Total output length will be multiplier * length(x)
        %   mode      : 'pre', 'post', or 'both' (default is 'both')
        %
        %   Example:
        %       x = [1 2 3];
        %       y = zero_pad_by_multiplier(x, 2, 'post'); % Pads to length 6
        
            if nargin < 3
                mode = 'both';
            end
        
            x = x(:);  % force column vector
            L = length(x);
            target_length = round(multiplier * L);
        
            if target_length < L
                error('Multiplier must be >= 1.');
            end
        
            pad_total = target_length - L;
        
            switch lower(mode)
                case 'pre'
                    x_padded = [zeros(pad_total, 1); x];
                case 'post'
                    x_padded = [x; zeros(pad_total, 1)];
                case 'both'
                    pre = floor(pad_total / 2);
                    post = ceil(pad_total / 2);
                    x_padded = [zeros(pre, 1); x; zeros(post, 1)];
                otherwise
                    error('Invalid mode. Use ''pre'', ''post'', or ''both''.');
            end
        
            % Restore row orientation if original was row vector
            if isrow(x)
                x_padded = x_padded.';
            end
        end



    end
end
