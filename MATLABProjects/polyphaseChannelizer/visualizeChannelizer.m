function visualizePolyphaseChannelizer(nChannels, nDec, fOffsetNorm, lenData, SNRdB, chanForSpectrogram)
% @brief visualization utility for the polyphaseChannelizer class. It mirrors the
%   unit-test stimulus, runs the DUT, and plots diagnostic views.
% @details Plots:
%   (1) Full-band input spectrum
%   (2) Per-channel output spectra (stacked) with expected tone offset marker
%   (3) Spectrogram of a selected output channel
%   (4) Prototype and polyphase subfilter responses
% @param nChannels: number of channels to produce from IQ (default 8)
% @param nDec: decimation ratio used by the channelizer (default 6)
% @param fOffsetNorm: fraction of channel bandwidth for the tone offset (default 0.1)
% @param lenData: input length (default 2^15)
% @param SNRdB: input SNR in dB (default 20)
% @param chanForSpectrogram: 1-based channel index for spectrogram (default 1)
% @return none
% @date 11/12/2025
% @author Daniel Beletsky

    %% ---------------------------- Defaults --------------------------------
    if nargin < 1 || isempty(nChannels),          nChannels = 8;      end
    if nargin < 2 || isempty(nDec),               nDec = 6;           end
    if nargin < 3 || isempty(fOffsetNorm),        fOffsetNorm = 0.1;  end
    if nargin < 4 || isempty(lenData),            lenData = 2^15;     end
    if nargin < 5 || isempty(SNRdB),              SNRdB = 20;         end
    if nargin < 6 || isempty(chanForSpectrogram), chanForSpectrogram = 1; end

    %% --------------------- Prototype Filter (matches tests) ----------------
    % Rough half-band per channel LP, equiripple with 1/f stopband shaping
    Fpass = 0.7*1/nChannels/2;
    Fstop = 1/nChannels/2;
    Apass = 0.5;
    Astop = 40;

    h  = fdesign.lowpass('Fp,Fst,Ap,Ast', Fpass, Fstop, Apass, Astop);
    Hd = design(h, 'equiripple', 'StopbandShape', '1/f', 'StopbandDecay', 3.25);
    proto = Hd.Numerator(:).';
    groupDelay = floor(length(proto)/2);

    %% ------------------------ DUT Construction ----------------------------
    % @note: uses your class as-is
    dut = polyphaseChannelizer(proto, groupDelay, nChannels, nDec);

    %% --------------------------- Stimulus ---------------------------------
    % Multi-tone: one per channel, each with same fractional offset within its BW
    n = 0:(lenData-1);
    chanCenters = (0:nChannels-1)/nChannels; % normalized to Fs
    tones = exp(1j*2*pi.*(chanCenters(:) + fOffsetNorm/(2*nChannels)).*n);
    x = sum(tones, 1);

    % Add white noise
    snrLin = db2mag(SNRdB);
    noise  = (randn(size(x)) + 1j*randn(size(x)))/(sqrt(2)*snrLin);
    x = x + noise;

    % Interleave I/Q as expected by DUT (I,Q,I,Q,...)
    xIQ = reshape([real(x); imag(x)], [], 1);

    %% --------------------------- Run DUT ----------------------------------
    Y = dut.filter(xIQ); % [time x nChannels], complex

    %% --------------------------- Analysis ---------------------------------
    % Input spectrum (context)
    NfftIn = 4*2^nextpow2(lenData);
    Xin = fftshift(fft(x, NfftIn));
    fIn = linspace(-0.5, 0.5, NfftIn);

    % Per-channel spectra (FFT along time)
    NfftCh = 2^nextpow2(size(Y,1));
    FY     = fft(Y, NfftCh, 1);
    FYmag  = abs(FY);
    fCh    = (0:(NfftCh-1))/NfftCh; % normalized 0..1 on decimated grid

    % Expected per-channel tone offset (0..1 within channel)
    % Derivation (same as in earlier helper): f_offset_ch = fOffsetNorm * nChannels/(2*nDec)
    expectedOffset = mod(fOffsetNorm * (nChannels/(2*nDec)), 1);
    expectedBin    = round(expectedOffset * NfftCh);

    %% ----------------------------- Plots ----------------------------------
    clf; shg;

    % (1) Full-band spectrum
    subplot(2,2,1);
    plot(fIn, 20*log10(abs(Xin)+1e-12), 'LineWidth', 1);
    grid on; axis tight;
    xlabel('Normalized frequency (full band)');
    ylabel('Magnitude (dB)');
    title(sprintf('Input Spectrum (N=%d)', lenData));
    hold on;
    for k = 0:nChannels
        xline((k/nChannels)-0.5, ':');
    end
    hold off;

    % (2) Per-channel spectra (stack first up-to-6 channels for readability)
    subplot(2,2,2);
    numToShow = min(nChannels, 6);
    hold on;
    for k = 1:numToShow
        pk = 20*log10(FYmag(:,k)+1e-12);
        plot(fCh + (k-1), pk, 'LineWidth', 1); % stack by offsetting x-axis per channel
        xExp = (expectedBin)/NfftCh + (k-1);
        plot([xExp xExp], [min(pk), max(pk)], '--');
    end
    hold off; grid on; axis tight;
    xlabel('Channel index (stacked)');
    ylabel('|FFT| (dB)');
    title(sprintf('Per-Channel Spectra (first %d); expected offset = %.3f', numToShow, expectedOffset));

    % (3) Spectrogram of selected channel
    subplot(2,2,3);
    ch = max(1, min(nChannels, chanForSpectrogram));
    win  = 1024;
    hop  = win/4;
    nfft = 2048;
    [S,F,T] = stft(double(Y(:,ch)), 'Window', hann(win,'periodic'), ...
        'OverlapLength', win-hop, 'FFTLength', nfft);
    imagesc(T, F, 20*log10(abs(S)+1e-12)); axis xy;
    xlabel('Time (frames)');
    ylabel('Normalized freq (0..1)');
    colorbar;
    title(sprintf('Channel %d Spectrogram', ch));

    % (4) Prototype + polyphase subfilter responses
    subplot(2,2,4);
    [Hproto, W] = freqz(proto, 1, 8192, 'whole'); W = W/(2*pi); % 0..1
    plot(W, 20*log10(abs(Hproto)+1e-12), 'LineWidth', 1); hold on;

    [fb, ~] = dut.getFilterBank(); % tapsPer x nChannels
    Np = 8192;
    for k = 1:min(nChannels,16) % limit visibility
        Hk = fft(fb(:,k), Np);
        plot((0:Np-1)/Np, 20*log10(abs(Hk)+1e-12), ':');
    end
    hold off; grid on; axis tight;
    xlabel('Normalized frequency (0..1)');
    ylabel('Magnitude (dB)');
    title(sprintf('Prototype & Polyphase Subfilters (showing %d)', min(nChannels,16)));
    legend({'Prototype','Subfilters (subset)'}, 'Location','southwest');

    %% ------------------------- Console Summary ----------------------------
    fprintf('\nVisualization summary:\n');
    fprintf('  nChannels = %d, nDec = %d, fOffsetNorm = %.3f, lenData = %d, SNR = %.1f dB\n', ...
        nChannels, nDec, fOffsetNorm, lenData, SNRdB);
    fprintf('  Expected per-channel tone offset (0..1) = %.6f (bin %d / %d)\n\n', ...
        expectedOffset, expectedBin, NfftCh);

end
