function channels = getElementData(M, SIG, methodIdx)
%getElementData Takes beamforming matrix and signal, outputs channel-specific
%data for coherence calculations.
%
%   Input
%   ~~~~~~~~~~~~~~~~~~
%   M         : Beamforming matrix [nSamples x length(SIG)]
%   SIG       : Signal column vector [nChannels*(nSamples+1) x 1]
%   methodIdx : Interpolation method index
%                 1: nearest, 2: linear (default), 3: quadratic,
%                 4: lanczos3, 5: 5points, 6: lanczos5
%
%   Output
%   ~~~~~~~~~~~~~~~~~~
%   channels  : Channel-specific data [nSamples x nChannels]

    %% Input validation
    if nargin < 3
        error('getElementData requires M, SIG, and methodIdx.');
    end
    if size(SIG, 2) > 1
        error('SIG must be a column vector. Check dimensions before calling getElementData.');
    end

    %% Precompute constants
    nSamples  = size(M, 1);
    nChannels = length(SIG) / nSamples;
    channels  = zeros(nSamples, nChannels);

    % Add small offset once outside loop to avoid zero-crossing suppression
    SIG = SIG + 1e-11;
    SIG_T = SIG';  % Transpose once; avoids repeated transposition in loop

    %% Main loop over samples
    for iSample = 1:nSamples
        % Extract non-zero weighted samples for this pixel
        vec = nonzeros(M(iSample, :) .* SIG_T);

        n = length(vec);

        % Accumulate interpolation method strides
        nOut  = floor(n / methodIdx);
        vec2  = vec(1:methodIdx:nOut*methodIdx);  % base stride

        for i = 2:methodIdx
            vec2 = vec2 + vec(i:methodIdx:nOut*methodIdx);
        end

        % Remove the offset contribution (methodIdx samples summed per output)
        vec2 = vec2 - 1e-11 * methodIdx;

        channels(iSample, 1:length(vec2)) = vec2;
    end

end