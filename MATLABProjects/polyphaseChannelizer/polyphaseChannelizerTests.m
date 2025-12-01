classdef polyphaseChannelizerTests < matlab.unittest.TestCase
    % @brief automated unit testing for the polyphase channelizer
    % @date 11/11/2025
    % @author Daniel Beletsky
    %

    %% Internal properties for passing between class setup and test
    properties
        polyphaseChan;
        filter;
    end

    %% Class Setup
    properties(ClassSetupParameter)
        nDec        = {2; 3; 4; 6; 7; 8; 14; 15; 16; 32; 64;};
        nChannels   = {4; 4; 4; 8; 8; 8; 16; 16; 16; 64; 64;};
    end

    methods(TestClassSetup, ParameterCombination = 'sequential')
        function classSetup(this, nDec, nChannels)
            % @brief setup the polyphase channelizer class. most of the 
            %   configuration is done in the constructor so we itterate through 
            %   those here.
            % @param nDec: decimation ratio to use in the channelizer
            % @param nChannels: number of channels to produce from IQ
            % @return none
            %

            % Define filter (just roughly passband)
            Fpass = 0.7*1/nChannels/2; % Passband Frequency
            Fstop = 1/nChannels/2; % Stopband Frequency
            Apass = 0.5;
            Astop = 40;
            h = fdesign.lowpass('Fp,Fst,Ap,Ast', Fpass, Fstop, Apass, Astop);
            Hd = design(h, 'equiripple', 'StopbandShape', '1/f', 'StopbandDecay', 3.25, 'SystemObject', true);
            this.filter = Hd.Numerator;
            groupDelay = floor(length(this.filter)/2);

            % Initialize channelizer
            this.polyphaseChan = polyphaseChannelizer(this.filter, groupDelay, nChannels, nDec);
        end
    end

    %% Method setup
    properties(MethodSetupParameter)
        seed = {0, 123, 98765};
    end

    methods(TestMethodSetup)
        function methodSetup(this, seed)
            % passed
        end
    end

    %% TestParameters and methods
    properties(TestParameter)
        fOffsetNorm = {0,0.1,0.3,-0.1,-0.3};
        lenData     = {1024, 5000, 2^17};
    end

    methods(Test)
        function testFilter(this, seed, nChannels, nDec, fOffsetNorm, lenData)
            % @brief test the filter function with various placed tones in
            %   each channel
            % @param seed: seed for rng in noise generation
            % @param nChannels: number of channels the channelizer makes
            % @param nDec: decimation ratio channelizer uses
            % @param fOffsetNorm: portion of channel bandwidth to place
            %   tone
            % @param lenData: length of data in
            % @return none
            %

            % Import relative comparisons 
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            % Make data
            chanCenters = [0:nChannels-1] / nChannels;
            tones = exp(1j*2*pi*[0:lenData-1].*(chanCenters(:) + fOffsetNorm/nChannels/2));
            dataIn = sum(tones, 1);
            rng(seed);
            snr = db2mag(20);
            noise = 1/snr*(randn(1,length(dataIn))/2 + 1j*randn(1,length(dataIn))/2);
            dataIn = dataIn + noise;

            % Filter Data
            dataInRe = real(dataIn);
            dataInIm = imag(dataIn);
            dataInterleaved = [dataInRe(:).'; dataInIm(:).'];
            dataInterleaved = dataInterleaved(:);
            dataOut = this.polyphaseChan.filter(dataInterleaved);

            % Check tones in right frequency bin
            [~,offset] = max(abs(fft(dataOut, [], 1)),[],1);
            offset = (offset-1)/size(dataOut,1);

            % Account for fftshift and decimation
            if any(offset > 0.5)
                offset = offset - 1;
            end
            offset = offset*2/nDec*nChannels;

            % Verify location of each tone
            for ii = 1:length(offset)
                this.verifyThat(offset(ii), IsEqualTo(fOffsetNorm, "Within", AbsoluteTolerance(2/size(dataOut,1))));
            end
        end

        function testGetters(this, nDec, nChannels)
            % @brief test all getting functions from the class
            % @param nDec: decmation ration channelizer is configured with
            % @param nChannels: number of channels channelizer produces
            % @return none
            %

            % Get the filter bank
            [fb, gd] = this.polyphaseChan.getFilterBank();

            % Verify filter bank sized right and group delay is right
            this.verifySize(fb, [ceil(length(this.filter)/nChannels), nChannels]);
            this.verifyEqual(gd, ceil(floor(length(this.filter)/2)/nDec));

            % Get the number of channels
            chan = this.polyphaseChan.getNChan();

            % Verify number channels correct
            this.verifyEqual(chan, nChannels);

        end
    end
end