%% One-Time Phase Calibration Script
% Author: Joe Pizzimenti

%% Run this ONCE with the Transmitter at where you want 0 degrees to be.
% It saves the phase offset to 'phase_cal.mat'.

clear; clc; close all;

% Config
addpath('C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\include\'); %AGAIN, THIS MAY BE DIFFERENT FOR DIFFERENT COMPUTERS!!!
uri = 'ip:192.168.2.1';
fc = 1.56e9; fs = 20e6;

% Reference Setup
rng(42); nSubcarriers = 64; cpLen = 16;
active_indices = [7:32 34:59];
known_freq = zeros(nSubcarriers, 1);
known_freq(active_indices) = sign(randn(length(active_indices), 1));
ref_time = ifft(ifftshift(known_freq)) * sqrt(nSubcarriers);
ref_time_cp = [ref_time(end-cpLen+1:end); ref_time];

% Setup Receiver
rx = adi.AD9361.Rx;
rx.uri = uri; rx.CenterFrequency = fc; rx.SamplingRate = fs;
rx.EnabledChannels = [1 2]; rx.kernelBuffersCount = 4;
rx.GainControlModeChannel0 = 'manual'; rx.GainChannel0 = 55;
rx.GainControlModeChannel1 = 'manual'; rx.GainChannel1 = 55;

disp('Collecting 50 packets for calibration...');
phase_diffs = [];

while length(phase_diffs) < 50
    data = rx();
    [xc, lags] = xcorr(data(:,1), ref_time_cp);
    [maxVal, maxIdx] = max(abs(xc));

    if maxVal < 1.5, continue; end % Skip noise

    startIdx = lags(maxIdx) - length(ref_time_cp) + 1;
    if startIdx < 1 || (startIdx + 80) > length(data), continue; end

    % Extract CSI
    symbol_start = startIdx + cpLen;
    rx_sym_time = data(symbol_start : symbol_start+nSubcarriers-1, :);
    rx_sym_freq = fftshift(fft(rx_sym_time));

    H1 = rx_sym_freq(active_indices, 1) ./ known_freq(active_indices);
    H2 = rx_sym_freq(active_indices, 2) ./ known_freq(active_indices);

    % Calculate Phase Difference between RX1 and RX2
    delta_phi = angle(mean(H2 ./ H1));
    phase_diffs = [phase_diffs; delta_phi];

    fprintf('.');
end

% Average the collected phase offsets
calibration_offset = mean(phase_diffs);
save('phase_cal.mat', 'calibration_offset');

disp(' ');
disp(['Calibration Complete! Offset: ' num2str(rad2deg(calibration_offset)) ' degrees']);
disp('You can now run the main Radar script.');
rx.release();
