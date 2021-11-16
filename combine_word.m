fname = {'beta1.mat' ; 'delta1.mat' ; 'theta1.mat'; 'alpha1.mat'};
combine_word = [];

for i = 1:4
    name = fname{i,:};
    eeg = load(sprintf('data/EEGdata/Synthetic EEG Fs512Hz/Noise/%s', name)).noisy_EEGsig;
    combine_word = [combine_word eeg];
end
