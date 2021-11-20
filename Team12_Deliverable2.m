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
    fname = 'beta2.mat';
    eeg = load(sprintf('data/EEGdata/Synthetic EEG Fs512Hz/Noise/%s', fname)).noisy_EEGsig;
else
    % Noiseless
    i = 1;
    eeg = load(sprintf('data/EEGdata/Synthetic EEG Fs512Hz/Noiseless/eeg%d.mat', i)).eeg;
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
plot(decision, 'LineWidth', 2)
title('Results');
xlabel('Frame Number');
ylabel('Decision');
legend('1 = delta, 2 = theta, 3 = alpha,  4 = beta');

%% C

% main{
% 
%     
%     decision(i+1) = make_decision(current_frame)
%     
%     store last 3 decisions
% 
%     decide what to play
% 
%     play word
% 
% 
% }

%
    % Pseudo code
%     
%     time_required = length(alphaWord)/fs_word; % duration of word (seconds)
%     frames_needed_to_play = time_required / hop;
%     tic
%     if (playing)
%             decision_check
%             save = frames_needed_to_play  + i
%     end
% 
%     % save is a future time step 
%     % once time step is reached, check past 3 decisions
%     % play word based on majority, 
%     % if all 3 = 0 || no distinct majority (etc), play no word
%     
%     clock
%     pause(computation)







