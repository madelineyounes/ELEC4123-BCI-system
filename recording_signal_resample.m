%% Resave as csv

alpha1out = table2array(alpha1out);
alpha2out = table2array(alpha2out);
beta1out = table2array(beta1out);
beta2out = table2array(beta2out);
theta1out = table2array(theta1out);
theta2out = table2array(theta2out);
delta1out = table2array(delta1out);
delta2out = table2array(delta2out);

csvwrite('recording_circuit_outputs/alpha1out.csv', alpha1out);
csvwrite('recording_circuit_outputs/alpha2out.csv', alpha2out);
csvwrite('recording_circuit_outputs/beta1out.csv', beta1out);
csvwrite('recording_circuit_outputs/beta2out.csv', beta2out);
csvwrite('recording_circuit_outputs/theta1out.csv', theta1out);
csvwrite('recording_circuit_outputs/theta2out.csv', theta2out);
csvwrite('recording_circuit_outputs/delta1out.csv', delta1out);
csvwrite('recording_circuit_outputs/delta2out.csv', delta2out);

%% Resample output from task 1 recording circuit to 512Hz

names = ["alpha1"; "alpha2"; "theta1"; "theta2"; "beta1"; "beta2"; "delta1"; "delta2"];

for i = 1:length(names)
%     x_ref = csvread('data/EEGdata/Synthetic EEG Fs512Hz/Noise/theta2.csv');

    data = csvread('recording_circuit_outputs/' + names(i) + 'out.csv');
    Tx = data(:, 1);
    x = data(:, 2);

    desiredFs = 512;
    [y, Ty] = resample(x, Tx, desiredFs);

%     figure(1);
%     subplot(2, 1, 1);
%     plot((0:length(x_ref)-1)/desiredFs, x_ref);
%     title("Raw Signal");
% 
%     subplot(2, 1, 2);
%     plot(Tx, x, '.-', Ty, y, 'o-');
%     legend('Original', 'Resampled');
%     title("Amplified Signal");

    csvwrite('recording_circuit_outputs/' + names(i) + 'out_512Hz.csv', y);
end

