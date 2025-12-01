Fs = 1;
N = 100;
f0 = 0.15;

t = (0:N-1)/Fs;
s = sin(f0 * 2 * pi * t);
w = blackman(N).';             % window
S = fft(s .* w);
S_mag = abs(S);
S_phase = angle(S);
f = (0:N-1)*Fs/N;

figure;
subplot(3,1,1); plot(t,s); title('Signal'); xlabel('Time (s)'); grid on;
subplot(3,1,2); plot(f,S_mag); title('Magnitude Spectrum'); xlabel('Frequency (Hz)'); grid on;
subplot(3,1,3); plot(f,S_phase); title('Phase Spectrum'); xlabel('Frequency (Hz)');

% One-sided spectrum (better for real signals)
P2 = abs(S)/N;
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f1 = (0:N/2)*(Fs/N);

figure;
subplot(3,1,1); plot(t, s); grid on; xlabel('Time (s)'); ylabel('s(t)'); title('Signal');
subplot(3,1,2); plot(f1, P1, '.-'); grid on; xlabel('Frequency (Hz)'); ylabel('|S(f)|'); title('Magnitude (one-sided)');
subplot(3,1,3); plot(f1, unwrap(angle(S(1:N/2+1))), '.-'); grid on; xlabel('Frequency (Hz)'); ylabel('Phase (rad)'); title('Phase (one-sided)');
