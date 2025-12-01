function show_ofdm_spectrum()
% show_ofdm_spectrum.m
% Displays FFT spectrum (frequency bins) of OFDM_TV_IQ_Image_Fs48khz.wav

wavFile = 'IQSig_96KHz.wav';
[y, fs] = audioread(wavFile);

% Use the same I/Q mapping (Q + jI)
I = y(:,1);
Q = y(:,2);
x = Q + 1i*I;

% Compute FFT of the whole file (or just a section if very long)
Nfft = 65536;
X = fftshift(fft(x, Nfft));
f = (-Nfft/2:Nfft/2-1)*(fs/Nfft);

figure('Name','OFDM Frequency Spectrum');
plot(f/1e3, 20*log10(abs(X)/max(abs(X)))); % normalize to dB
grid on;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB, normalized)');
title('OFDM Spectrum (Frequency Bins)');
xlim([-fs/2 fs/2]/1e3);

% Optional: zoom in on the active OFDM band (~±6 kHz)
hold on; xline([-6 6],'--r'); hold off;

fprintf('Displayed %d FFT bins, fs = %.0f Hz.\n', Nfft, fs);
end
