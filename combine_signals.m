%% Combine noisy EEG signals

noisy_combined = [];
names = ["beta"; "delta"; "theta"; "alpha"];
% Corresponds to "I", "Love", "Electrical", "Engineering"
set = 1;

for i = 1:length(names)
    name = names(i) + num2str(set) + ".mat";
    % Transpose to become column vector
    eeg = load("data/EEGdata/Synthetic EEG Fs512Hz/Noise/" + name).noisy_EEGsig';
    noisy_combined = [noisy_combined; eeg];
end

csvwrite("data/EEGdata/Synthetic EEG Fs512Hz/Noise/noisy_combined" + num2str(set) + ".csv", noisy_combined);

%% Combine noiseless EEG signals

noiseless_combined = [];

for i = 1:6
    % Transpose to become column vector
    eeg = load("data/EEGdata/Synthetic EEG Fs512Hz/Noiseless/eeg" + num2str(i) + ".mat").eeg';
    noiseless_combined = [noiseless_combined; eeg];
end

csvwrite("data/EEGdata/Synthetic EEG Fs512Hz/Noiseless/noiseless_combined.csv", noiseless_combined);
