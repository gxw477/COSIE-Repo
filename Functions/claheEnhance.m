function enhanced = claheEnhance(img, varargin)
% CLAHEENHANCE  Contrast-Limited Adaptive Histogram Equalization (CLAHE)
%               for B-mode ultrasound images.
%
%   enhanced = claheEnhance(img)
%   enhanced = claheEnhance(img, Name, Value, ...)
%
% -------------------------------------------------------------------------
%   INPUTS (required)
% -------------------------------------------------------------------------
%   img        - Input B-mode image.
%                Accepted classes: uint8, uint16, double, single.
%                Grayscale only (2-D array).
%
% -------------------------------------------------------------------------
%   NAME-VALUE PARAMETERS (all optional)
% -------------------------------------------------------------------------
%   'NumTiles'    - 1x2 vector [rows, cols] specifying the number of
%                   contextual tiles the image is divided into.
%                   Each tile computes its own histogram.
%                   Default: [8, 8].
%                   Guideline: larger tiles → more global; smaller → more
%                   local. Tile size should be at least 3–4x the expected
%                   feature diameter.
%
%   'ClipLimit'   - Normalised clip limit in (0, 1].
%                   Controls contrast amplification ceiling.
%                   Excess histogram counts are redistributed uniformly.
%                   Default: 0.02.
%                   Lower → less contrast, less noise amplification.
%                   Higher → stronger contrast, more noise risk.
%
%   'NBins'       - Number of histogram bins per tile.
%                   Default: 256.
%
%   'Distribution'- Desired histogram shape after equalization.
%                   'uniform'     : flat histogram (classic HE). [default]
%                   'rayleigh'    : right-skewed; often better for US
%                                   speckle (brighter mean).
%                   'exponential' : aggressive dark-region boost.
%
%   'Alpha'       - Shape parameter for 'rayleigh' and 'exponential'
%                   distributions. Default: 0.4.
%                     rayleigh    : alpha controls Rayleigh sigma
%                     exponential : alpha controls decay rate
%
%   'Range'       - Output intensity range.
%                   'original' : map to [min(img), max(img)]. [default]
%                   'full'     : map to full class range
%                                (e.g. [0 255] for uint8).
%
% -------------------------------------------------------------------------
%   OUTPUT
% -------------------------------------------------------------------------
%   enhanced   - Enhanced image, same class and size as input.
%
% -------------------------------------------------------------------------
%   ALGORITHM OVERVIEW
% -------------------------------------------------------------------------
%   1. Normalise input to [0, 1].
%   2. Pad image so it divides evenly into NumTiles.
%   3. For each tile:
%        a. Compute histogram with NBins bins.
%        b. Clip histogram at ClipLimit * (tile_pixels / NBins) and
%           redistribute clipped counts uniformly across all bins.
%        c. Compute the CDF and map it through the chosen Distribution
%           to form a tile transfer function.
%   4. Bilinear interpolation between the four nearest tile transfer
%      functions to produce a smooth, artifact-free result (no tile
%      boundary seams).
%   5. Crop padding, cast to original class.
%
% -------------------------------------------------------------------------
%   EXAMPLES
% -------------------------------------------------------------------------
%   % Default — good general-purpose starting point
%   out = claheEnhance(img);
%
%   % Finer tiles + lower clip: conservative, low-noise enhancement
%   out = claheEnhance(img, 'NumTiles', [16 16], 'ClipLimit', 0.01);
%
%   % Rayleigh distribution: better preserves ultrasound speckle texture
%   out = claheEnhance(img, 'Distribution', 'rayleigh', 'Alpha', 0.3);
%
%   % Strong local contrast for deep tissue boundary visualisation
%   out = claheEnhance(img, 'NumTiles', [4 4], 'ClipLimit', 0.05, ...
%                      'Distribution', 'exponential');
%
% -------------------------------------------------------------------------
%   REFERENCES
% -------------------------------------------------------------------------
%   Zuiderveld, K. (1994). "Contrast Limited Adaptive Histogram
%   Equalization." Graphics Gems IV, Academic Press, pp. 474–485.
%
%   Pizer, S.M. et al. (1987). "Adaptive histogram equalization and its
%   variations." Computer Vision, Graphics, and Image Processing, 39(3).
%
% -------------------------------------------------------------------------
%   See also: adapthisteq, histeq, sigmoidEnhance.

    % ================================================================== %
    %  Parse inputs
    % ================================================================== %
    p = inputParser();
    p.FunctionName = mfilename;

    addRequired(p,  'img');
    addParameter(p, 'NumTiles',     [8 8],      @(x) isvector(x) && numel(x)==2 && all(x>=2));
    addParameter(p, 'ClipLimit',    0.02,       @(x) isscalar(x) && x>0 && x<=1);
    addParameter(p, 'NBins',        256,        @(x) isscalar(x) && x>=2 && x==floor(x));
    addParameter(p, 'Distribution', 'uniform',  @(x) ismember(x,{'uniform','rayleigh','exponential'}));
    addParameter(p, 'Alpha',        0.4,        @(x) isscalar(x) && x>0);
    addParameter(p, 'Range',        'original', @(x) ismember(x,{'original','full'}));

    parse(p, img, varargin{:});
    opts = p.Results;

    numTiles  = double(opts.NumTiles(:)');   % [nTR, nTC]
    clipLimit = opts.ClipLimit;
    nBins     = opts.NBins;
    distrib   = opts.Distribution;
    alpha     = opts.Alpha;
    rangeMode = opts.Range;

    % ================================================================== %
    %  Validate image
    % ================================================================== %
    validateattributes(img, {'uint8','uint16','double','single'}, ...
        {'nonempty','2d','real','nonsparse'}, mfilename, 'img', 1);

    origClass = class(img);

    % ================================================================== %
    %  Normalise to [0, 1]
    % ================================================================== %
    switch origClass
        case 'uint8'
            x = double(img) / 255;
            fullMin = 0; fullMax = 255;
        case 'uint16'
            x = double(img) / 65535;
            fullMin = 0; fullMax = 65535;
        otherwise   % double / single
            x = double(img);
            if max(x(:)) > 1 || min(x(:)) < 0
                x = (x - min(x(:))) / (max(x(:)) - min(x(:)) + eps);
            end
            fullMin = 0; fullMax = 1;
    end

    origMin = min(x(:));
    origMax = max(x(:));

    [nR, nC] = size(x);
    nTR = numTiles(1);   % tile rows
    nTC = numTiles(2);   % tile cols

    % ================================================================== %
    %  Pad image so it divides evenly into tiles
    % ================================================================== %
    tileH = ceil(nR / nTR);   % tile height in pixels
    tileW = ceil(nC / nTC);   % tile width  in pixels

    padR = tileH * nTR - nR;
    padC = tileW * nTC - nC;

    % Reflect-pad to avoid edge artefacts
    xPad = padarray(x, [padR, padC], 'symmetric', 'post');
    [nRp, nCp] = size(xPad);

    % ================================================================== %
    %  Pre-compute tile transfer functions
    %  tileTF(:, r, c) = transfer function for tile (r,c)
    % ================================================================== %
    tileTF = zeros(nBins, nTR, nTC);

    binEdges  = linspace(0, 1, nBins + 1);   % nBins+1 edges → nBins bins
    clipCount = clipLimit * (tileH * tileW);  % absolute clip threshold

    for tr = 1:nTR
        for tc = 1:nTC
            % Extract tile
            rowIdx = (tr-1)*tileH + (1:tileH);
            colIdx = (tc-1)*tileW + (1:tileW);
            tile   = xPad(rowIdx, colIdx);

            % Histogram
            counts = histcounts(tile(:), binEdges);   % 1 x nBins

            % Clip & redistribute
            counts = clahe_clip(counts, clipCount);

            % CDF → normalised to [0,1]
            cdf = cumsum(counts);
            cdf = (cdf - cdf(1)) / (cdf(end) - cdf(1) + eps);

            % Map CDF through desired distribution
            tileTF(:, tr, tc) = cdf2dist(cdf, distrib, alpha, nBins);
        end
    end

    % ================================================================== %
    %  Bilinear interpolation of transfer functions across tiles
    % ================================================================== %
    enhanced_pad = zeros(nRp, nCp);

    % Centre coordinates of each tile (in padded image pixel coords)
    tileCentreR = ((0:nTR-1) + 0.5) * tileH;   % 1 x nTR
    tileCentreC = ((0:nTC-1) + 0.5) * tileW;   % 1 x nTC

    for r = 1:nRp
        for c = 1:nCp
            pix = xPad(r, c);

            % Bin index for this pixel
            binIdx = max(1, min(nBins, ceil(pix * nBins)));
            if binIdx == 0; binIdx = 1; end

            % Find surrounding tile centres
            [t1r, t2r, wr] = interp_weights(r, tileCentreR, nTR);
            [t1c, t2c, wc] = interp_weights(c, tileCentreC, nTC);

            % Bilinear blend of the four surrounding tile TFs
            v  = (1-wr)*(1-wc) * tileTF(binIdx, t1r, t1c) ...
               + (1-wr)*   wc  * tileTF(binIdx, t1r, t2c) ...
               +    wr *(1-wc) * tileTF(binIdx, t2r, t1c) ...
               +    wr *   wc  * tileTF(binIdx, t2r, t2c);

            enhanced_pad(r, c) = v;
        end
    end

    % ================================================================== %
    %  Crop padding and map to output range
    % ================================================================== %
    enh = enhanced_pad(1:nR, 1:nC);   % remove padding

    % Map [0,1] to desired output range
    switch rangeMode
        case 'original'
            enh = enh * (origMax - origMin) + origMin;
        case 'full'
            % already in [0,1]; scale to full class range below
    end

    % Cast back to original class
    switch origClass
        case 'uint8'
            if strcmp(rangeMode,'full')
                enhanced = uint8(enh * 255);
            else
                enhanced = uint8(enh * 255);
            end
        case 'uint16'
            if strcmp(rangeMode,'full')
                enhanced = uint16(enh * 65535);
            else
                enhanced = uint16(enh * 65535);
            end
        otherwise
            enhanced = cast(enh, origClass);
    end
end


% ======================================================================= %
%  LOCAL FUNCTION: clip histogram and redistribute excess counts
% ======================================================================= %
function counts = clahe_clip(counts, clipCount)
% Clip each bin at clipCount and redistribute the excess uniformly.
    excess      = sum(max(0, counts - clipCount));
    counts      = min(counts, clipCount);
    redistribute = floor(excess / numel(counts));
    remainder    = mod(excess, numel(counts));

    counts = counts + redistribute;

    % Distribute remainder to first R bins
    if remainder > 0
        counts(1:remainder) = counts(1:remainder) + 1;
    end
end


% ======================================================================= %
%  LOCAL FUNCTION: map CDF to desired output distribution
% ======================================================================= %
function tf = cdf2dist(cdf, distrib, alpha, nBins)
% Transform a normalised CDF [0,1] into a transfer function according to
% the desired output histogram distribution.
    switch distrib
        case 'uniform'
            % Classic HE: TF = CDF directly
            tf = cdf;

        case 'rayleigh'
            % Inverse Rayleigh CDF: x = sigma * sqrt(-2*ln(1-p))
            % Clamp to avoid log(0)
            p  = min(cdf, 1 - eps);
            tf = alpha * sqrt(-2 * log(1 - p + eps));
            tf = tf / max(tf(:) + eps);   % normalise to [0,1]

        case 'exponential'
            % Inverse Exponential CDF: x = -ln(1-p) / lambda
            p  = min(cdf, 1 - eps);
            tf = -log(1 - p + eps) / alpha;
            tf = tf / max(tf(:) + eps);   % normalise to [0,1]
    end

    tf = max(0, min(1, tf));   % hard clamp
end


% ======================================================================= %
%  LOCAL FUNCTION: find two nearest tile centres and interpolation weight
% ======================================================================= %
function [i1, i2, w] = interp_weights(pos, centres, n)
% Given a pixel position `pos` and a vector of tile centre coordinates
% `centres`, returns the indices of the two nearest centres (i1, i2)
% and the weight w toward i2 (so weight toward i1 is 1-w).

    if pos <= centres(1)
        i1 = 1; i2 = 1; w = 0;
    elseif pos >= centres(end)
        i1 = n; i2 = n; w = 0;
    else
        i2 = find(centres >= pos, 1, 'first');
        i1 = i2 - 1;
        w  = (pos - centres(i1)) / (centres(i2) - centres(i1) + eps);
    end
end