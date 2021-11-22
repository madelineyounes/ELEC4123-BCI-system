%% Noiseless to wav

clc; clear;

name = 'eeg';
fs = 512;

for i = 1:6
    eeg = load([name num2str(i) '.mat']).eeg;
    audiowrite([name num2str(i) '.wav'],eeg, fs);  
end

%% Noiseless to csv

clc; clear;

for i = 1:6
    % Transpose to become column vector
    eeg = load(['eeg' num2str(i) '.mat']).eeg';
    csvwrite(['eeg' num2str(i) '.csv'], eeg);  
end

%% Noisy to csv

clc; clear;

name = {'alpha1';
        'alpha2';
        'theta1';
        'theta2';
        'beta1';
        'beta2';
        'delta1';
        'delta2'};

for i = 1:length(name)
    % Transpose to become column vector
    eeg = load([name{i} '.mat']).noisy_EEGsig';
    csvwrite([name{i} '.csv'], eeg);
end

%% Real to csv

clear;

data = load('data/EEGdata/Real EEG data/rec1.mat');
fs = data.fs;
eeg = data.channelData(:, 1);   % Get data in the first channel (column)
csvwrite('data/EEGdata/Real EEG data/rec1.csv', eeg);


