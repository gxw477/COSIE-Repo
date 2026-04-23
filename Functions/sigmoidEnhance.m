function enhanced = sigmoidEnhance(img, gain, cutoff)
% SIGMOIDENHANCE  Sigmoid contrast enhancement for B-mode ultrasound images.
%
%   enhanced = sigmoidEnhance(img)
%   enhanced = sigmoidEnhance(img, gain, cutoff)
%
%   INPUTS:
%     img     - Input B-mode image (uint8, uint16, or double [0,1]).
%     gain    - Steepness of the sigmoid curve (default: 10).
%               Higher values → sharper contrast transition.
%     cutoff  - Normalised intensity threshold at which sigmoid = 0.5
%               (default: 0.5). Range (0, 1).
%               Lower values → brighten darker tissues;
%               Higher values → suppress more background noise.
%
%   OUTPUT:
%     enhanced - Enhanced image, same class as input.
%
%   DESCRIPTION:
%     Applies a pixel-wise sigmoid transfer function:
%
%         S(x) = 1 / (1 + exp(-gain * (x - cutoff)))
%
%     then linearly rescales the output to span [0, 1] so that the
%     full dynamic range is preserved after enhancement.
%
%   EXAMPLES:
%     % Basic usage with defaults
%     img = dicomread('bmode.dcm');
%     out = sigmoidEnhance(img);
%     imshow(out);
%
%     % Aggressive contrast stretch (high gain, low cutoff)
%     out = sigmoidEnhance(img, 15, 0.35);
%
%     % Gentle enhancement preserving speckle texture
%     out = sigmoidEnhance(img, 6, 0.5);
%
%   NOTE:
%     For a stack of frames (cine-loop), loop over the third dimension:
%       for k = 1:size(volume,3)
%           volume(:,:,k) = sigmoidEnhance(volume(:,:,k), gain, cutoff);
%       end
%
%   See also: imadjust, adapthisteq, histeq.

    % ------------------------------------------------------------------ %
    %  Input validation & defaults
    % ------------------------------------------------------------------ %
    if nargin < 2 || isempty(gain);   gain   = 10;  end
    if nargin < 3 || isempty(cutoff); cutoff = 0.5; end

    validateattributes(gain,   {'numeric'}, {'scalar','real','positive'},          mfilename, 'gain',   2);
    validateattributes(cutoff, {'numeric'}, {'scalar','real','>',0,'<',1},         mfilename, 'cutoff', 3);
    validateattributes(img,    {'uint8','uint16','double','single'}, {'nonempty'}, mfilename, 'img',    1);

    % ------------------------------------------------------------------ %
    %  Normalise to [0, 1]
    % ------------------------------------------------------------------ %
    origClass = class(img);

    switch origClass
        case 'uint8'
            x = double(img) / 255;
        case 'uint16'
            x = double(img) / 65535;
        otherwise                          % double / single — assume [0,1]
            x = double(img);
            if max(x(:)) > 1              % handle arbitrary double range
                x = (x - min(x(:))) / (max(x(:)) - min(x(:)) + eps);
            end
    end

    % ------------------------------------------------------------------ %
    %  Sigmoid transfer function
    % ------------------------------------------------------------------ %
    % S(x) = 1 / (1 + exp(-gain*(x - cutoff)))
    S = 1 ./ (1 + exp(-gain .* (x - cutoff)));

    % ------------------------------------------------------------------ %
    %  Contrast-stretch output to full [0, 1] dynamic range
    % ------------------------------------------------------------------ %
    Smin = min(S(:));
    Smax = max(S(:));
    S    = (S - Smin) / (Smax - Smin + eps);

    % ------------------------------------------------------------------ %
    %  Cast back to original data type
    % ------------------------------------------------------------------ %
    switch origClass
        case 'uint8'
            enhanced = uint8(S * 255);
        case 'uint16'
            enhanced = uint16(S * 65535);
        otherwise
            enhanced = cast(S, origClass);
    end
end