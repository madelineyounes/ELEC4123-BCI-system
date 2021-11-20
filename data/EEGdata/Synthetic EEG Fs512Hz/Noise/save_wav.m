%% Most recently used to convert noiseless eeg to .wav

clc; clear;

name = 'eeg';
fs = 512;

for i = 1:6
    eeg = load([name num2str(i) '.mat']).eeg;
    audiowrite([name num2str(i) '.wav'],eeg, fs);  
end
