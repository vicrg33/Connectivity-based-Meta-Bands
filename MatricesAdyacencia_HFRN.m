clc;
clear;
close all;
ft_defaults;

%% PARAMETERS
signal_type = 'MEG';
force_params = 0; % 1 if you want to calculate all the parameters of each script on the same iteration (e.g. if PLI is selected, PLI, wPLI and PLV will be calculated)
trial_overlapping_or = 0; % Trial overlapping (in percentage)
trial_len_or = 5; % Trial length (in seconds)
max_trials = []; % Maximum number of trials. Empty = use all trials
fs_vec = []; % Sampling frequencies to be tested. Empty = do not modify the sampling frequency
mode = 2; % 0 = no filtering; 1 = bands; 2 = window
    % Filtering parameters
    filter_orders = [500]; % Filter order
    bands = 'classic'; % For mode 1: 'classic' (classical bands) or 'custom' (user-defined bands)
    band_init = [1 15 30]; % For mode 1: user defined bands (initial cutoff frequencies)
    bands_end = [15 30 70]; % End cutoff frequencies
    freqs_overlapping_or = [50]; % For mode 2: frequency overlapping values (in percentage) to be tested
    freqs_resolution = [1]; % For mode 2: frequency resolution, i.e. frequency window size values (in Hz) to be tested
    f_start = 1; % For mode 0 or 2: lower frequency limit (in Hz)
    f_end = 70; % For mode 0 or 2: higher frequency limit (in Hz)
        
    
%% PATH CONFIGURATION

% Paths with files
fs = % Original fs
files = dir('Path/To/Files.mat');

% If no sampling frequencies are selected, use the default of the database
if isempty(fs_vec)
    fs_vec = fs;
end 

%% LOOP HEADERS    

for fs_new = fs_vec % Loop through all the selected sampling frequencies...
    for filter_order = filter_orders
        for freq_resolution = freqs_resolution % Loop through all the selected frequential resolution...
            for freq_overlapping_or = freqs_overlapping_or % Loop through all the selected frequency overlapping values...
%% CONFIGURATION

                % Configure the parameters and define the bands
                % Filter design
                Den = 1; % Den
                flag = 'scale';  % Sampling Flag
                win = hamming(filter_order+1); % Window      

                % Bands design                 
                freq_overlapping = 1 - freq_overlapping_or/100;
                bands_init = f_start:freq_resolution*freq_overlapping:(f_end-freq_resolution);                                       
                bands_end = bands_init + freq_resolution;                  
                freq_size = length(bands_init);

                % Other parameters
                trial_overlapping = 1 - trial_overlapping_or/100;
                trial_len = trial_len_or * fs_new; % Trial length fron sec to samples


                %% PROCESSING
                % For all the subjects
                for n_file  = length(files):-1:1
                    
                    load([files(n_file).folder, '/', files(n_file).name]); % Load the data file (signal)
                    cleansignals.data = signal;

                    % Re-referencing (Only for EEG signals)
                    if isequal(signal_type,'EEG')
                        signal_average = []; matriz_ref = [];
                        signal_average = nanmean(cleansignals.data,2);
                        matriz_ref = cleansignals.data - repmat(signal_average,1,size(cleansignals.data,2));
                        cleansignals.data = [];
                        cleansignals.data = matriz_ref;
                    end                    

                    % Subsample to the desired sampling frequency
                    cleansignals.data = resample(cleansignals.data,fs_new,fs);                    

                    % Normalize signals
                    cleansignals.data = zscore(cleansignals.data);                    

                    % More parameters...
                    num_chann = size(cleansignals.data,2); % Channel or ROI number
                    trial_size = trial_len:floor(trial_len*trial_overlapping):size(cleansignals.data,1);
                    trial_size = trial_size - trial_len + 1;
                    % Remove N first and last trials to avoid filter
                    % transient effects
                    n_trial_discard = ceil((filter_order/2) / trial_len); % Filter has a transient of filter_order/2 (actually less, due to reflection of filtfilt)
                    trial_size = trial_size(1+n_trial_discard:end-n_trial_discard);

                    % Only "max_trials" trials
                    if ~isempty(max_trials)
                        trial_limit = min(max_trials,length(trial_size));
                        if (trial_limit < max_trials) % If the subject do not have enough trials, skip it
                            continue
                        end
                    else
                        trial_limit = length(trial_size); 
                    end
                    trial_size = trial_size(1:trial_limit);

                    adj_matrix = NaN(freq_size,num_chann,num_chann,length(trial_size));
                    
                    for n_freq = 1:freq_size % Loop through all the frequency bins

                        % Filter if a filter is needed
                        % Creating the filter for the current frequency bin
                        Fc1  = bands_init(n_freq); % First cutoff freq
                        Fc2  = bands_end(n_freq); % Second cutoff freq
                        Num  = fir1(filter_order, [Fc1 Fc2]/(fs_new/2), 'bandpass', win, flag); % Bandpass Design
                        % Filter the signal
                        signal_filtered = filtfilt(Num, Den, cleansignals.data);

                        for n_trial = 1:length(trial_size) % Loop through all the trials

                            % Select only the current trial
                            signal_trial = signal_filtered(trial_size(n_trial):trial_size(n_trial)+trial_len-1, :);

                            % Ortogonalizate the signal (For AEC_ort calculation)
                            signal_ort = ortogonalization_optimizada(signal_trial);

                            parfor n_chann = 1:num_chann % For each channel
                                % Of note, is not the same ortogonalize channel A
                                % regarding B that channel B regarding A. Thus, the
                                % matrix will not be symmetrical. To address this issue
                                % is common to average both connectivity measures (A>B
                                % and B>A)

                                % Calcualte AEC
                                AEC_tmp = calculo_AEC(squeeze(signal_ort(:, :, n_chann)));
                                adj_matrix(n_freq, :, n_chann, n_trial) = AEC_tmp(:, n_chann); % Store it in the adjacecy matrix (This values will be replaced)
                            end

                            % Averega both connectivity measures
                            AEC1 = triu(squeeze(adj_matrix(n_freq,:,:,n_trial)));
                            AEC2 = tril(squeeze(adj_matrix(n_freq,:,:,n_trial))).';
                            AECglobal = (AEC1+AEC2)/2;
                            adj_matrix(n_freq,:,:,n_trial) = abs(triu(AECglobal, 1)+AECglobal.'); % Store it in the adjacency matrix

                        end
                    end

                end
            end
        end
    end
end
