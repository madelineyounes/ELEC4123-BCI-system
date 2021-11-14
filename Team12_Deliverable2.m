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

clear; close all;
fs = 512;

% Load EEG signal
noisy = true;
if noisy
    % Noisy
    fname = 'alpha1.mat';
    eeg = load(sprintf('data/EEGdata/Synthetic EEG Fs512Hz/Noise/%s', fname)).noisy_EEGsig;
else
    % Noiseless
    i = 2;
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

%% Plot Truncated EGG signal
truncate1 = 11 * 0.125 * fs; % Truncate signal
truncate2 = 22 * 0.125 * fs;
eeg = eeg(truncate1:truncate2);
L = length(eeg);

figure(3);
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

%% Signal processing
% Derived from DP DSP Task 4 starting code

% Load EEG signal
noisy = true;
if noisy
    % Noisy
%     fname = 'beta2.mat';
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
win_offset = 0.8;
win_scale = 0.2;
frame_win = kaiser(N)';
delta_win = (win_scale*hamming(5) + win_offset)';
theta_win = (win_scale*hamming(5) + win_offset)';
alpha_win = (win_scale*hamming(5) + win_offset)';
beta_win = (win_scale*hamming(29) + win_offset)';

% Noise Parameters
noise_mean = 0;
noise_std = 0;

% Noise Learning Rate
mean_eta = 0.9; 
std_eta = 0.9;

% Loop through all frames in the signal
for i = 1:num_frames
    j = (i-1)*hop + 1;
    frame = eeg(j:j+N-1);
    
    % Frame windowing
%     frame = frame .* frame_win;
    
    % FFT (check magnitude vs. power spectrum)
    Frame = abs(fft(frame)).^2;
    Frame = Frame(1:(N/2)); % Single sided magnitude response
    
    % Frequency band windowing
    Frame_delta = Frame(1:5);   % 0 - 4 Hz 
    Frame_theta = Frame(5:9);   % 4 - 8 Hz
    Frame_alpha = Frame(9:13);  % 8 - 12 Hz
    Frame_beta = Frame(13:41);  % 12 - 40 Hz
    
    Frame_delta = Frame_delta .* delta_win;
    Frame_theta = Frame_theta .* theta_win;
    Frame_alpha = Frame_alpha .* alpha_win;
    Frame_beta = Frame_beta .* beta_win;
    
    % Adaptive threshold
    Frame_noise = Frame(42:256);
    if i == 1
        noise_mean = mean(Frame_noise);
        noise_std = std(Frame_noise);
    else
        noise_mean = noise_mean * mean_eta + mean(Frame_noise)*(1-mean_eta);
        noise_std = noise_std * std_eta + std(Frame_noise)*(1-std_eta);
    end
    
    noise_threshold = noise_mean + 2*noise_std;
    delta_avg = mean(Frame_delta) - noise_threshold;
    theta_avg = mean(Frame_theta) - noise_threshold;
    alpha_avg = mean(Frame_alpha) - noise_threshold;
    beta_avg = mean(Frame_beta) - noise_threshold;
    
    % Estimate noise from 40 - 255 Hz
    % Update noise estimate every frame
   
    if (all([delta_avg, theta_avg, alpha_avg, beta_avg] <= 0))
        decision(i) = 0;
    else 
        [~,index] = max([delta_avg theta_avg alpha_avg beta_avg]);
        decision(i) = index;
    end 
    
%     figure(4)
%     plot(Frame);
%     yline(delta_avg, 'r', 'LineWidth', 1.5);
%     yline(theta_avg, 'g', 'LineWidth', 1.5);
%     yline(alpha_avg, 'b', 'LineWidth', 1.5);
%     yline(beta_avg, 'y', 'LineWidth', 1.5);
%     yline(noise_threshold, 'LineWidth', 1.5);
%     legend('Power Spectrum' ,'delta' ,'theta', 'alpha', 'beta', 'noise threshold');
%     xlim([0 45]);
%     pause(2)
end

% xline(4, 'r', 'LineWidth', 1.5);    % Delta waves: 0–4 Hz


figure(5)
plot(decision)
title('Results');
xlabel('Frame Number');
ylabel('Decision');
legend('1 = delta, 2 = theta, 3 = alpha,  4 = beta');





