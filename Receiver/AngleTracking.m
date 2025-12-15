%% SpotFi Radar Test (Atempts to Ignore Multipath)
% Author: Joe Pizzimenti

%% WARNING! MUST RUN Calibration.m before running this for calibration prior to angle tracking
% Uses Weighted Linear Fit to attempt to ignore multipath
% A sloped phase vs subcarrier plot means mulipath is too strong for angle
% tracking with 2 antennas

clear; clc; close all;

% Load Calibration file output from Calibration.m
if exist('phase_cal.mat', 'file')
    load('phase_cal.mat');
    disp(['Loaded Calibration Offset: ' num2str(rad2deg(calibration_offset)) ' deg']);
else
    warning('No calibration file found! Using 0.');
    calibration_offset = 0;
end

% Config
addpath('C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\include\'); %THIS WILL BE DIFFERENT DEPENDING ON HOW LIBIIO IS INSTALLED!!!
uri = 'ip:192.168.2.1';
fc = 1.56e9; fs = 20e6; c = 3e8;
antSpacing = c/fc/2; %about 0.0962

scanAngles = -90:1:90;

% Setup same RNG seed as transmitter (reference setup)
rng(42); nSubcarriers = 64; cpLen = 16;
active_indices = [7:32 34:59];
known_freq = zeros(nSubcarriers, 1);
known_freq(active_indices) = sign(randn(length(active_indices), 1));
ref_time = ifft(ifftshift(known_freq)) * sqrt(nSubcarriers);
ref_time_cp = [ref_time(end-cpLen+1:end); ref_time];

% Setup Receiver
rx = adi.AD9361.Rx; %CONFIGURED AS A CUSTOM AD9361 (not a standard pluto chip)
rx.uri = uri; rx.CenterFrequency = fc; rx.SamplingRate = fs;
rx.EnabledChannels = [1 2]; rx.kernelBuffersCount = 4;
rx.GainControlModeChannel0 = 'manual'; rx.GainChannel0 = 55;
rx.GainControlModeChannel1 = 'manual'; rx.GainChannel1 = 55;

disp('Starting Radar...');
figure('Name', 'Angle Tracker', 'Position', [100 100 900 700]);

% Create plots
subplot(2,1,1);
hPlot = plot(scanAngles, zeros(size(scanAngles)), 'LineWidth', 2);
grid on; xlim([-90 90]); ylim([-30 0]);
xlabel('Angle (deg)'); ylabel('Spectrum (dB)');
title('Angle Estimate');
xline_h = xline(0, '--r', 'Est Target');

subplot(2,1,2);
hPhase = plot(nan, nan, '.'); hold on;
hFit = plot(nan, nan, 'r', 'LineWidth', 2);
grid on; title('Phase Difference vs Subcarrier (Slope = Multipath)');
xlabel('Subcarrier Index'); ylabel('Phase Diff (rad)');
legend('Raw Diff', 'Robust Fit');

while true
    data = rx();
    [xc, lags] = xcorr(data(:,1), ref_time_cp);
    [maxVal, maxIdx] = max(abs(xc));

    if maxVal < 1.5, continue; end

    startIdx = lags(maxIdx) - length(ref_time_cp) + 1;
    if startIdx < 1 || (startIdx + 80) > length(data), continue; end

    % Extract CSI
    symbol_start = startIdx + cpLen;
    rx_sym_time = data(symbol_start : symbol_start+nSubcarriers-1, :);
    rx_sym_freq = fftshift(fft(rx_sym_time));

    H1 = rx_sym_freq(active_indices, 1) ./ known_freq(active_indices);
    H2 = rx_sym_freq(active_indices, 2) ./ known_freq(active_indices);

    % Phase Estimation
    raw_diff = H2 ./ H1;
    raw_diff = raw_diff * exp(-1j * calibration_offset);
    weights = abs(H1) .* abs(H2);
    weights = weights.^2; % Square it to punish weak signals heavily (weighted squares)
    phases = angle(raw_diff);
    phases = unwrap(phases); % Remove 2*pi jumps

    robust_phase = sum(phases .* weights) / sum(weights);

    % Update Debug Plot (Bottom) (SLOPE MEANS MULTIPATH IS PRESENT!!!)
    set(hPhase, 'XData', 1:length(phases), 'YData', phases);
    set(hFit, 'XData', [1 length(phases)], 'YData', [robust_phase robust_phase]);
    H_clean = [1; exp(1j * robust_phase)]; % Reconstructed 2x1 vector

    % Scan Angles
    spectrum = zeros(1, length(scanAngles));
    for i = 1:length(scanAngles)
        a = [1; exp(-1j * 2 * pi * antSpacing * sind(scanAngles(i)) * (fc/c))];
        spectrum(i) = abs(a' * H_clean)^2; % Beamforming Power (Bartlett is more stable here than MUSIC for 2 antennas)
    end

    spectrum = 10*log10(spectrum / max(spectrum)); % normalize

    % Visualization
    [~, peakIdx] = max(spectrum);
    est_angle = scanAngles(peakIdx);

    set(hPlot, 'YData', spectrum);
    set(xline_h, 'Value', est_angle);
    title(hPlot.Parent, sprintf('Target: %d deg', est_angle));

    drawnow limitrate;
end
