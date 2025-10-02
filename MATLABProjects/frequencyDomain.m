Fs = 1;
N = 100;

t = (0:N-1)/Fs;
s = sin(0.15 * 2 * pi * t);
S = fft(s);
S_mag = abs(S);
S_phase = angle(S);
f = (0:N-1)*Fs/N;

figure;
subplot(3,1,1); plot(t,s); title('Signal'); xlabel('Time (s)'); grid on;
subplot(3,1,2); plot(f,S_mag); title('Magnitude Spectrum'); xlabel('Frequency (Hz)'); grid on;
subplot(3,1,3); plot(f,S_phase); title('Phase Spectrum'); xlabel('Frequency (Hz)');
