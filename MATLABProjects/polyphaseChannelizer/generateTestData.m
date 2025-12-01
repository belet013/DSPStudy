function testData = generateTestData(seed, nDec, nChan)
% @brief simple function generating test input/output for simulink,
%    need to make it just for matlab.
% @param seed: Random seed
% @param kGain: Gain to use for the filter
%
% @return testData: structure of test data
%

% Verify the inputs are there or provide default values
    if ~exist('seed','var') || isempty(seed)
        seed = 1234567;
    end
     if ~exist('nDec','var') || isempty(nDec)
        seed = 4;
    end
     if ~exist('seed','nChan') || isempty(nChan)
        seed = 4;
    end

    % Make a random stream