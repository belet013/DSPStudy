function decode_ofdm_tv()
% Final OFDM-TV decoder tuned to OFDM_TV_IQ_Image_Fs48khz.wav
% Key details:
%   fs = 48 kHz
%   Chirp: 100 ms, BW = 12 kHz, positive slope
%   Correct I/Q packing: x = Q + j*I
%   CRITICAL reshape fix: row-major emulation in MATLAB

% ----------------------------
% 0) Load WAV
% ----------------------------
wavFile = 'OFDM_TV_IQ_Image_Fs48khz.wav';
fprintf('Loading: %s\n', wavFile);
[y, fs] = audioread(wavFile);
if abs(fs - 48e3) > 1
    warning('Unexpected fs = %.2f Hz (expected 48 kHz). Proceeding...', fs);
end
if size(y,2) < 2
    error('Expected stereo I/Q WAV. Found %d channel(s).', size(y,2));
end

% For THIS file, channels are swapped: use Q + j*I
I = y(:,1);
Q = y(:,2);
x = Q + 1i*I;

% ----------------------------
% 1) Build chirp & matched filter
% ----------------------------
dt    = 1/fs;
pw    = 100e-3;      % 100 ms
bw    = 12e3;
slope = bw/pw;
t     = (dt:dt:pw).';
t     = t - pw/2;
lfm   = exp(1i*pi*slope.*t.^2);    % positive-slope chirp
h     = conj(flipud(lfm(:)));

% Use I channel for timing (works robustly)
r = conv(I, h, 'full');
[~, peakIdx] = max(abs(r));
L = length(h);
startIdx = peakIdx - L + 1;               % first sample AFTER chirp (1-based)
startIdx = max(1, startIdx);

% ----------------------------
% 2) Slice the OFDM block
% ----------------------------
numRows = 480;
numCols = 1024;
Nneeded = numRows * numCols;

if startIdx + Nneeded - 1 > length(x)
    % Clamp to last valid start if needed
    startIdx = length(x) - Nneeded + 1;
    if startIdx < 1
        error('File too short: need %d complex samples; have %d.', Nneeded, length(x));
    end
end

ofdm = x(startIdx : startIdx + Nneeded - 1);

% ----------------------------
% 3) Reshape (ROW-MAJOR emulation!), FFT, phase-diff image
% ----------------------------
% NumPy row-major equivalent in MATLAB:
Xtime = reshape(ofdm, numCols, numRows).';   % <--- THE IMPORTANT FIX

Xf = fftshift(fft(Xtime, [], 2), 2);         % FFT per row, then shift

% Phase difference between adjacent rows
pic = angle( Xf(2:end,:) ./ Xf(1:end-1,:) ); % (479 x 1024)

% ----------------------------
% 4) Plots (MF sanity + images)
% ----------------------------
figure('Name','Matched filter magnitude');
plot(abs(r)); grid on;
xlabel('Sample index (conv domain)'); ylabel('|r|');
title('Matched filter output magnitude (I channel)');
hold on; xline(peakIdx,'--r','Peak'); hold off;

figure('Name','Phase-difference image (as-is)');
imagesc(pic); axis image; colormap gray; colorbar;
title('Phase-difference image (as-is)');

figure('Name','Phase-difference image (negated)');
imagesc(-pic); axis image; colormap gray; colorbar;
title('Phase-difference image (negated)');

fprintf('Done. Payload start index: %d\n', startIdx);
end
