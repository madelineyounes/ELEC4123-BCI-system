%% View real EEG data

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
    fname = 'theta2.mat';
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

%% Plot Truncated EEG Signal
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
