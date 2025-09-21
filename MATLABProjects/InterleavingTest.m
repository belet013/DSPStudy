clear; close all; clc;

Fs = 48e3; % sample rate
T = 0.5; % duration (s)
N = round(T*Fs); % samples
fA = 1e3; % stream A frequency
fB = 100; % stream B frequency
ampA = 1; % stream A amplitude
ampB = 0.1; % stream B amplitude

Nfft = 1024; % FFT length (vector size)
hop = Nfft/4; % hop size (75% overlap)
window = blackman(Nfft); %window for FFT frames
 % making two streams (scalars per time sample)
t = (0:N-1)/Fs;
streamA = ampA*cos(2*pi*fA*t); % stream A
streamB = ampB*cos(2*pi*fA*t); % stream B
 %create a vector of the two streams. (Parallel samples per time sample)
vec2 = [streamA; streamB];
 % Plot the two input streams.
figure;
subplot(3,1,1); plot(t,streamA); title('Stream A'); grid on; xlim([0 5e-3]);
subplot(3,1,2); plot(t,streamB); title('Stream B'); grid on; xlim([0 5e-3]);
subplot(3,1,3); plot(t, vec2(1,:),'b'); hold on; plot(t,vec2(2,:), 'r'); grid on; xlim([0 5e-3]); legend('A','B');
title('2-element Vector at each time (stacked for display)'); 
 %Interleave data like [A(1) B(1) A(2) B(2)]
interleaved_data = reshape(vec2, 1, []); % size 1 x 2N;
 % Plot interleaved data vs original
figure;
subplot(3,1,1); stem(0:15, streamA(1:16),'filled'); title('Stream A (first 16 samples)'); grid on; 
subplot(3,1,2); stem(0:15, streamB(1:16),'filled'); title('Stream B (first 16 samples)'); grid on;
subplot(3,1,3); stem(0:31, interleaved_data(1:32),'filled'); title('Interleaved data, first 32?'); grid on;
 % Stream to vectors for block processing in fft. Frame overlapping with % buffer, columns are frames of length Nfft buffer(x, Nfft, Noverlap, % 'nodelay')
Noverlap = Nfft-hop; % overlap
framed = buffer(interleaved_data, Nfft, Noverlap, 'nodelay'); % Nfft x nFrames
nFrames = size(framed,2);
 % Apply window and compute FFT (vector operation per frame)
X = fft(framed.*window, Nfft, 1);
mag = 20*log10(abs(X(1:Nfft/2+1,:))+1e-12);
 % Gives us spectogram
f = (0:Nfft/2)*Fs/Nfft;
tau = (0:nFrames-1)*hop/Fs;
figure;
imagesc(tau,f/1e3,mag); axis xy; colorbar; xlabel('Time(s)'); ylabel('Freq(kHz)');
 % De-interleave or serialize back into 2 streams
vec2_recovered = reshape(interleaved_data,2,[]);
A_rec = vec2_recovered(1,:);
B_rec = vec2_recovered(2,:);
 % Verify it worked
fprintf('Max |A-A_rec| = %.3g\n', max(abs(streamA-A_rec)));
fprintf('Max |B-B_rec| = %.3g\n', max(abs(streamB-B_rec)));

figure;
subplot(2,1,1); plot(t, A_rec); title('Recovered Stream A'); grid on; xlim([0 5e-3]);
subplot(2,1,2); plot(t, B_rec); title('Recovered Stream B'); grid on; xlim([0 5e-3]);