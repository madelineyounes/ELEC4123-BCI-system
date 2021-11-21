%% Function to decide which word to play
% Inputs: EEG frame, noise_mean, noise_std
% Output: EEG Class (delta/theta/alpha/beta)

% Assume window parameters are global
function [decision, noise_mean, noise_std] = make_decision(frame, noise_mean, noise_std )
    
    % Frame windowing
%     frame = frame .* frame_win;

    % Noise Learning Rate
    mean_eta = 0.2; 
    std_eta = 0.2;

    % Windows (row vector)
    win_offset = 0.8;
    win_scale = 0.2;
    delta_win = (win_scale*hamming(5) + win_offset)';
    theta_win = (win_scale*hamming(5) + win_offset)';
    alpha_win = (win_scale*hamming(5) + win_offset)';
    beta_win = (win_scale*hamming(29) + win_offset)';

    % FFT (Obtain power spectrum)
    Frame = abs(fft(frame)).^2;
    Frame = Frame * (1/512);
    Frame = Frame(1:(length(frame)/2)); % Single sided power spectrum
    
    % Frequency band windowing
    Frame_delta = Frame(1:5);   % 0 - 4 Hz 
    Frame_theta = Frame(5:9);   % 4 - 8 Hz
    Frame_alpha = Frame(9:13);  % 8 - 12 Hz
    Frame_beta = Frame(13:41);  % 12 - 40 Hz
    
    Frame_delta = Frame_delta .* delta_win;
    Frame_theta = Frame_theta .* theta_win;
    Frame_alpha = Frame_alpha .* alpha_win;
    Frame_beta = Frame_beta .* beta_win;
    
    % Estimate noise from 40 - 255 Hz
    Frame_noise = Frame(42:256);
    if (noise_mean == 0 && noise_std == 0)
        noise_mean = mean(Frame_noise);
        noise_std = std(Frame_noise);
    else
        % Update noise estimate every frame
        noise_mean = noise_mean * mean_eta + mean(Frame_noise)*(1-mean_eta);
        noise_std = noise_std * std_eta + std(Frame_noise)*(1-std_eta);
    end
    
    noise_threshold = noise_mean + 1.4*noise_std;
    delta_avg = mean(Frame_delta) - noise_threshold;
    theta_avg = mean(Frame_theta) - noise_threshold;
    alpha_avg = mean(Frame_alpha) - noise_threshold;
    beta_avg = mean(Frame_beta) - noise_threshold;
    
    if (all([delta_avg, theta_avg, alpha_avg, beta_avg] <= 0))
        decision = 0;
    else
        [~, index] = max([delta_avg theta_avg alpha_avg beta_avg]);
        decision = index;
    end 

     % Iteratively updating figure
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

   