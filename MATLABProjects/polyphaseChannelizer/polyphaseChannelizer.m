classdef polyphaseChannelizer < handle %#codevisualizeChannelizegen
    % @brief polyphase channelizer code for taking in IQ symbols and producing a
    %       configurable number of p channels.
    % @note a good resource for this is: 
    %       Polyphase Channelizer Demystified
    % @note a good resource on non-maximally decimated implementations is:
    %       https://www.mathworks.com/help/dsp/ref/channelizer.html
    %   see resoruce [2] in this link. This paper describes the implementation
    %   in greater detail
    % @date: 11/10/2025
    % @author: Daniel Beletsky
    %

    %% Public Member Variables
    properties(Access = public)
    end

    %% Private Member Variables
    properties(Access = private)
        nChan;
        filterBank;
        tapsPerFilter;
        ppGroupDelay;
        nDec;
        maxDec              = false;
        nSizeSection        = uint32(0);
        nSectionsIn         = uint32(0);
        nStates             = uint32(0);
        oversampleFactor    = single(0);
        stateStepSize       = uint32(0);
    end

    %% Class Setup Methods
    methods(Access = public)
        function this = polyphaseChannelizer(filter, groupDelay, nChan, nDec)
            % @brief: construtor for the polyphase channelizer. Takes a full rate 
            %   filter and its associated group delay and forms a polyphase 
            %   channelizer with nChan channels that decimates the data by nDec.
            % @param filter: array of tap values at full rate
            % @param nChan: # of channels to produce from a full rate stream
            % @param nDec: decimation rate to use. If nDec is equal to nChan then 
            %   the polyphase channelizer is maximally decimated.
            % @return polyphaseChannelizer object
            %

            % Ensure inputs are valid
            if nDec > nChan
                error("polyphaseChannelizer:inavlidInputs", strcat('Error \n The Decimation ration cannot exceed ', 'the number of channels.'));
            end

            % Determine Filter Bank dimensions
            this.nChan              = nChan;
            this.nDec               = nDec;
            nTaps                   = length(filter);
            this.tapsPerFilter      = ceil(nTaps/nChan);
            % Allocate Filter Bank
            this.filterBank         = zeros(this.tapsPerFilter,nChan);

            % Assign Filter Bank polyphase components
            for iChan = 1:this.nChan
                ppFilt = filter(iChan:nChan:end);
                len_ppFilt = length(ppFilt);
                this.filterBank(:,iChan) = [ppFilt(:); zeros(this.tapsPerFilter - len_ppFilt,1)];
            end
            this.ppGroupDelay = ceil(groupDelay/this.nDec); % group delay at output

            % Determine if we can use a maximally decimated architecture.
            this.maxDec = nChan == nDec;

            % Determine block sizes if not maximally decimated.
            if ~this.maxDec
                this.nSizeSection       = uint32(1);
                this.nSectionsIn        = uint32(this.nDec/this.nSizeSection);
                this.nStates            = uint32(lcm(this.nDec, this.nChan)/this.nDec);
                this.oversampleFactor   = single(this.nChan)/single(this.nDec);
                this.stateStepSize      = uint32((this.nChan - this.nDec)/gcd(this.nDec, this.nChan));
            end
        end
    end

    %% Setter/Getter Methods
    methods(Access = public)
        function [filterBank, groupDelays] = getFilterBank(this)
            % @brief get the filter bank and associated group delay in the taps
            % @return filterbank the nChan subfilters
            % @return groupDelay: group delay associated with each subfilter
            filterBank = this.filterBank;
            groupDelays = this.ppGroupDelay;
        end

        function nChans = getNChan(this)
            % @brief get the number of channels
            % @return nChans: number of channels to produceout of an IQ stream.
            nChans = this.nChan;
        end
    end

    %% Public Processing Methods  
    methods(Access = public)
            function iqChanOut = filter(this, iqFullBand)
            % @brief filter the fill band iq presented as interleaved I,Q data into
            %   individual channels with polyphase channelizer implementation
            % @param iqFullBand: Full band IQ data interleaved as I, Q, I, Q,....
            % @return iqChnOut: Decimated channels of size iTime x nChan
            %

            % Deinterleave the IQ data
            I_data = iqFullBand(1:2:end);
            Q_data = iqFullBand(2:2:end);
            iq = I_data(:) + 1j*Q_data(:);

            % Allocate Buffer
            lenData = length(iq);
            lenStream = ceil(lenData/this.nDec);
            dataStreams = complex(zeros(lenStream, this.nChan));

            if this.maxDec  
                % Do max decimated proccessing flow (no memory swaps)

                % Resample data streams
                for iChan = 1:this.nChan
                    ppComponent                             = iq(iChan:this.nChan:end);
                    len_ppComponent                         = length(ppComponent);
                    dataStreams(:, this.nChan - iChan + 1)  = [ppComponent; 
                        complex(zeros(lenStream - len_ppComponent, 1))];
                end

                % Filter    
                filteredLen = lenStream + this.tapsPerFilter - 1;
                filtData = complex(zeros(filteredLen, this.nChan));

                for iChan = 1:this.nChan
                    % Maybe make ths fast convolution
                    filtData(:, iChan) = conv(dataStreams(:, iChan),...
                        this.filterBank(:, iChan));
                end
                % Allocate output Buffer
                iqChanOut = complex(zeros(filteredLen - this.ppGroupDelay,...
                    this.nChan));
                
                % FFT Across time to unmix channels
                for iSamp = this.ppGroupDelay + 1:filteredLen
                    iqChanOut(iSamp - this.ppGroupDelay,:) = fft(fliplr(filtData(iSamp, :)));
                end

            else    
                % Compute indicies (could possibly change this to be on initialization)
                old_rows        = (1:this.nChan).'; % Fix this, makes VSCode look weird, take away ""
                new_rows        = mod(uint32(old_rows) - 1 + this.nSectionsIn, this.nChan) + 1; % +/- for Matlab indexing
                colShift        = uint32(old_rows > new_rows);
                allOldCols      = repmat(uint32((1:this.tapsPerFilter)), this.nChan,1);
                allNewCols      = allOldCols + colShift;
                allNewCols      = allNewCols(:);
                allOldCols      = allOldCols(:);
                allNewRows      = repmat(new_rows, this.tapsPerFilter,1);
                allOldRows      = repmat(old_rows, this.tapsPerFilter,1);
                logIndsToRemove = allNewCols > this.tapsPerFilter; % remove off end assignments
                allNewCols(logIndsToRemove) = [];       
                allOldCols(logIndsToRemove) = [];
                allNewRows(logIndsToRemove) = [];
                allOldRows(logIndsToRemove) = [];
                newInds         = sub2ind([this.nChan, this.tapsPerFilter], allNewRows, allNewCols);
                oldInds         = sub2ind([this.nChan, this.tapsPerFilter], allOldRows, allOldCols);

                % Do data alignment plus filtering (coupled)
                filteredLen = lenStream + uint32(ceil(single(this.tapsPerFilter)*this.oversampleFactor)) - 1;
                filtData    = complex(zeros(filteredLen, this.nChan));
                dataBuf     = complex(zeros(this.nChan, this.tapsPerFilter));
                nSampsIn    = this.nSizeSection * this.nSectionsIn;

                % Step through each time step 
                for iTime = 1:filteredLen
                    % Serpentine data through Buffer
                    dataBuf(newInds)    = dataBuf(oldInds);

                    % Move in data
                    endMoveIn   = iTime*nSampsIn;
                    startMoveIn = (iTime-1)*nSampsIn+1;
                    if startMoveIn > lenData
                        dataBuf(1:nSampsIn,1) = zeros(nSampsIn, 1);
                    elseif endMoveIn > lenData
                        dataBuf(1:nSampsIn,1) = [flipud(iq(startMoveIn:lenData));zeros(endMoveIn-lenData,1)];
                    else
                        dataBuf(1:nSampsIn,1) = flipud(iq((iTime-1)*nSampsIn + 1: iTime*nSampsIn));
                    end

                    % Form interproducts
                    filtData(iTime,:) = dot(this.filterBank.', dataBuf,2); % note not commutative
                end

                % Allocate output Buffer
                iqChanOut = complex(zeros(ceil(filteredLen - this.ppGroupDelay), this.nChan));

                % Unmix channel streams
                for iSamp = this.ppGroupDelay:filteredLen
                    % Cyclic shift to align polyphase
                    dataMat = reshape(filtData(iSamp,:).', [], this.nStates);
                    dataShifted = circshift(dataMat, int32(mod(this.stateStepSize*iSamp, this.nStates)),2);

                    % FFT
                    iqChanOut(iSamp-this.ppGroupDelay+1,:) = fft(flipud(dataShifted(:))).';
                end
            end
        end
    end

    %% Private Processing Methods
    methods(Access = private)
        function [realOut,imagOut] = cMultIn32(real1, imag1, real2, imag2)
            realOut = real1.*real2 - imag1.*imag2;
            imagOut = real1.*imag2 + imag1.*real2;
        end
    end
end



