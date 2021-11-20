fs = 512;

channelData = load('data/EEGdata/Synthetic EEG Fs512Hz/Noise/theta1.mat').noisy_EEGsig;
L = length(channelData);
time = (0:L-1)*(1/fs);
M = [time' channelData'];
writematrix(M, 'theta1.csv');