fname = {'beta1.mat' . 'delta1.mat' ; 'theta1.mat'; 'alpha1.mat'};


combine_word = [];

for i = 1:length(fname)
    [alphaWord, fs_word]= audioread("data/Word wav files/alphaWord.wav");
    name = fname{i,:};
    eeg = load(sprintf('data/EEGdata/Synthetic EEG Fs512Hz/Noise/%s', name)).noisy_EEGsig;
    combine_word = [combine_word eeg];
end

