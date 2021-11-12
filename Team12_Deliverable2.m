% ELEC4123 Elective - Analog Electronics & DSP
% Team 12

% Delta waves: 0–4 Hz
% Theta waves: 4–8 Hz
% Alpha waves: 8–12 Hz
% Beta waves: 12–40 Hz

%% View real EEG data (not on Github due to large file sizes)

clear;
i = 1;
load(sprintf('EEGdata/Real EEG data/rec%d.mat', i));
figure(1);
time = (0:length(channelData)-1)*(1/fs);
plot(time, channelData);
if exist('channelNames', 'var')
    legend(channelNames);
end
xlabel('Time (sec)');
ylabel('Amplitude (uV)');
grid on;

%% View synthetic EEG data

clear;
fs = 512;
noisy = true;
% Load EEG signal
if noisy
    % Noisy
    fname = 'alpha1.mat';
    eeg = load(sprintf('EEGdata/Synthetic EEG Fs512Hz/Noise/%s', fname)).noisy_EEGsig;
else
    % Noiseless
    i = 1;
    eeg = load(sprintf('EEGdata/Synthetic EEG Fs512Hz/Noiseless/eeg%d.mat', i)).eeg;
end
L = length(eeg);

% Plot EEG signal in time domain
figure(2);
subplot(2, 1, 1);
time = (0:L-1)*(1/fs);
plot(time, eeg*1e6);    % Assuming original data is in units of V
xlabel('Time (sec)');
ylabel('Amplitude (uV)');
grid on;

% Plot EEG signal in frequency domain
subplot(2, 1, 2);
freq = linspace(-fs/2, fs/2, L);
plot(freq, abs(fftshift(fft(eeg))));
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

