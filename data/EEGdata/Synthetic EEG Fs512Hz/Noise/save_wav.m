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
    eeg = load(['eeg' num2str(i) '.mat']).eeg;
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
    eeg = load([name{i} '.mat']).noisy_EEGsig';
    csvwrite([name{i} '.csv'], eeg);
end
