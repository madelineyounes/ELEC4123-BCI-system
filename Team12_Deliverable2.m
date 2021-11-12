% ELEC4123 Elective - Analog Electronics & DSP
% Team 12

%% View real EEG data (not on Github due to large file sizes)

clear;
i = 2;
load(sprintf('data/EEGdata/Real EEG data/rec%d.mat', i));
L = length(channelData);
figure(1);
sgtitle(sprintf('data/Real EEG Signal - rec%d.mat', i));

% Plot EEG signal in time domain
subplot(2, 1, 1);
time = (0:L-1)*(1/fs);
plot(time, channelData);
if exist('channelNames', 'var')
    legend(channelNames);
end
title('Time Domain Plot');
xlabel('Time (sec)');
ylabel('Amplitude (uV)');
grid on;

% Plot EEG signal in frequency domain
subplot(2, 1, 2);
freq = linspace(-fs/2, fs/2, L);
plot(freq, abs(fftshift(fft(channelData./1e6))));   % Convert from uV to V
if exist('channelNames', 'var')
    legend(channelNames);
end
title('Magnitude Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Annotate figure
xlim([0 45]);
xline(4, 'r', 'LineWidth', 1.5);    % Delta waves: 0–4 Hz
xline(8, 'r', 'LineWidth', 1.5);    % Theta waves: 4–8 Hz
xline(12, 'r', 'LineWidth', 1.5);   % Alpha waves: 8–12 Hz
xline(40, 'r', 'LineWidth', 1.5);   % Beta waves: 12–40 Hz
text(2, 0, 'Delta', 'color', 'r', 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(6, 0, 'Theta', 'color', 'r', 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(10, 0, 'Alpha', 'color', 'r', 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(26, 0, 'Beta', 'color', 'r', 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

%% View synthetic EEG data

clear;
fs = 512;

% Load EEG signal
noisy = false;
if noisy
    % Noisy
    fname = 'beta1.mat';
    eeg = load(sprintf('data/EEGdata/Synthetic EEG Fs512Hz/Noise/%s', fname)).noisy_EEGsig;
else
    % Noiseless
    i = 1;
    eeg = load(sprintf('data/EEGdata/Synthetic EEG Fs512Hz/Noiseless/eeg%d.mat', i)).eeg;
end
L = length(eeg);

% Plot EEG signal in time domain
figure(2);
if noisy
    sgtitle(sprintf('Noisy - %s', fname));
else
    sgtitle(sprintf('Noiseless - eeg%d.mat', i));
end 

subplot(2, 1, 1);
time = (0:L-1)*(1/fs);
plot(time, eeg*1e6);    % Assuming original data is in units of V
title('Time Domain Plot');
xlabel('Time (sec)');
ylabel('Amplitude (uV)');
grid on;

% Plot EEG signal in frequency domain
subplot(2, 1, 2);
freq = linspace(-fs/2, fs/2, L);
plot(freq, abs(fftshift(fft(eeg))));
title('Magnitude Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Annotate figure
xlim([0 45]);
xline(4, 'r', 'LineWidth', 1.5);    % Delta waves: 0–4 Hz
xline(8, 'r', 'LineWidth', 1.5);    % Theta waves: 4–8 Hz
xline(12, 'r', 'LineWidth', 1.5);   % Alpha waves: 8–12 Hz
xline(40, 'r', 'LineWidth', 1.5);   % Beta waves: 12–40 Hz
text(2, 0, 'Delta', 'color', 'r', 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(6, 0, 'Theta', 'color', 'r', 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(10, 0, 'Alpha', 'color', 'r', 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(26, 0, 'Beta', 'color', 'r', 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

%% Signal processing
% Derived from DP DSP Task 4 starting code

% Load EEG signal
noisy = false;
if noisy
    % Noisy
    fname = 'alpha1.mat';
    eeg = load(sprintf('data/EEGdata/Synthetic EEG Fs512Hz/Noise/%s', fname)).noisy_EEGsig;
else
    % Noiseless
    i = 1;
    eeg = load(sprintf('data/EEGdata/Synthetic EEG Fs512Hz/Noiseless/eeg%d.mat', i)).eeg;
end
L = length(eeg);
fs = 512;
N = 1 * fs;
hop = 0.125 * fs;
num_frames = ceil((L-N)/hop);
decision = zeros(num_frames, 1);

% Windows
frame_win = hamming(N);
frame_delta = hamming(4);
frame_theta = hamming(4);
frame_alpha = hamming(4);
frame_beta = hamming(28);

% Loop through all frames in the signal
for i = 1:num_frames
    j = (i-1)*hop + 1;
    frame = eeg(j:j+N-1);
    
    % Frame windowing
    frame = frame .* frame_win;
    
    % FFT (check magnitude vs. power spectrum)
    Frame = abs(fft(frame)).^2;
    Frame = Frame(1:(N/2)); % Single sided magnitude response
    
    % Frequency band windowing
    

    % Adaptive threshold
    
end





