%% Rescale audio files to be more friendly to the DAC

rescaled_max = 3;   % volts

[alphaWord, fs] = audioread('alphaWord.wav');
alphaWord = alphaWord * 3 / max(abs(alphaWord));
audiowrite('alphaWord_rescaled.wav', alphaWord, fs);

% Note: audiowrite() clips the signal outside the range [-1, 1]

% betaWord = betaWord * 3 / max(abs(betaWord));
% deltaWord = deltaWord * 3 / max(abs(deltaWord));
% thetaWord = thetaWord * 3 / max(abs(thetaWord));
