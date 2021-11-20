% ELEC4123 Elective - Analog Electronics & DSP
% Team 12

%% Signal processing
% Derived from DP DSP Task 4 starting code

% Load in words
[alphaWord, fs_word]= audioread("data/Word wav files/alphaWord.wav");
[betaWord, fs_word]= audioread("data/Word wav files/betaWord.wav");
[deltaWord, fs_word]= audioread("data/Word wav files/deltaWord.wav");
[thetaWord, fs_word]= audioread("data/Word wav files/thetaWord.wav");

% Load EEG signal
noisy = true;
if noisy
    % Noisy
    name = "beta2";
    eeg = load("data/EEGdata/Synthetic EEG Fs512Hz/Noise/" + name + ".mat").noisy_EEGsig;
else
    % Noiseless
    i = 1;
    eeg = load(sprintf("data/EEGdata/Synthetic EEG Fs512Hz/Noiseless/eeg%d.mat", i)).eeg;
end

% eeg = combine_word;

L = length(eeg);
fs = 512;
N = 1 * fs;
hop = 0.125 * fs; 
num_frames = ceil((L-N)/hop);
decision = zeros(num_frames+1, 1);
time_finish = datetime('now');

% Noise Parameters
noise_mean = 0;
noise_std = 0;

for i = (2:num_frames + 1) 
    tic
    
    % Frame choose
    j = (i-2)*hop + 1;
    frame = eeg(j:j+N-1);

    % Classify EEG into alpha/theta/delta/beta
    [decision(i), noise_mean, noise_std] = make_decision(frame, noise_mean, noise_std);
    
    % Play word
    if (datetime('now') >= time_finish)
        if (decision(i) == decision(i-1))
            if (decision(i) == 1)
                soundsc(deltaWord, fs_word);
                time_finish = datetime('now') + seconds(length(deltaWord)/fs_word);
            elseif (decision(i) == 2)
                soundsc(thetaWord, fs_word);
                time_finish = datetime('now') + seconds(length(thetaWord)/fs_word);
            elseif (decision(i) == 3)
                soundsc(alphaWord, fs_word);
                time_finish = datetime('now') + seconds(length(alphaWord)/fs_word);
            elseif (decision(i) == 4)
                soundsc(betaWord, fs_word);
                time_finish = datetime('now') + seconds(length(betaWord)/fs_word);
            end
        end
    end
    
    pause(hop/fs - toc)
end

figure(5)
plot(decision(2:end), 'LineWidth', 2)
title('Results');
xlabel('Frame Number');
ylabel('Decision');
legend('1 = delta, 2 = theta, 3 = alpha,  4 = beta');


%% Plot and play output from C++ code

% Load test results
% name = "theta1";
test_results_decisions = csvread("data/Test Results/" + name + "_decisions.csv");
test_results_audio_samples = csvread("data/Test Results/" + name + "_audio_samples.csv");

figure(6);
sgtitle("Test Results: " + name + ".csv");

subplot(2, 1, 1);
plot(test_results_decisions, 'LineWidth', 2);
title('Decisions');
xlabel('Frame Number');
ylabel('Decision');
legend('1 = delta, 2 = theta, 3 = alpha,  4 = beta');

subplot(2, 1, 2);
plot((0:(length(test_results_audio_samples)-1))/fs_word, test_results_audio_samples);
title('Audio Samples');
xlabel('Time (sec)');
ylabel('Amplitude');

% Play sound
soundsc(test_results_audio_samples, fs_word);
