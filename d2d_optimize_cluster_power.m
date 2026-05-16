function Hc = d2d_optimize_cluster_power(cl, i, grids, quant, cutoff, PG, cfg)
%
% This function determines the optimal scaling of a single cluster tensor
% such that its total power (after thresholding) matches the target cluster
% power as closely as possible.
%
% The optimization is performed by searching over a set of candidate
% dynamic range (DR) levels. For each candidate DR:
%   - A cluster tensor is generated
%   - A global threshold (cutoff) is applied
%   - The tensor is scaled accordingly
%   - The resulting power is compared with the target power
%
% The DR value that minimizes the power mismatch is selected.
%
% This step ensures that:
%   - the reconstructed cluster matches the statistical power model
%   - the truncation effects (due to DR) are properly accounted for
%
% ========================================================================
% INPUT
% ========================================================================
% cl     : cluster structure
%          - pwr   : cluster power (dB)
%          - DS, delay, DoA, DoD, etc.
%
% i      : cluster index
%
% grids  : precomputed grids (delay, DoA, DoD)
%
% quant  : quantization lookup table
%
% cutoff : global threshold (linear scale)
%          determined by dominant cluster and DR
%
% PG     : path gain (linear scale)
%          used to convert relative cluster power to absolute scale
%
% cfg    : configuration structure
%          - max_DR : maximum DR search range (in dB)
%
% ========================================================================
% OUTPUT
% ========================================================================
% Hc : scaled cluster tensor [N_delay × N_DoA × N_DoD]
%      best match to target power under DR constraint
%
% ========================================================================


% ------------------------------------------------------------------------
% TARGET POWER (LINEAR SCALE)
% ------------------------------------------------------------------------
% Convert cluster power (dB) + path gain to linear scale
target = 10^((cl.pwr(i) + 10*log10(PG)) / 10);


% ------------------------------------------------------------------------
% INITIALIZATION
% ------------------------------------------------------------------------
best = inf;    % best power mismatch so far
Hc = [];       % output tensor (empty if no valid candidate found)


% ------------------------------------------------------------------------
% SEARCH OVER DYNAMIC RANGE VALUES
% ------------------------------------------------------------------------
% Iterate over candidate DR levels (1 dB resolution)
for jj = 1:cfg.max_DR

    % ------------------------------------------------------------
    % Step 1: Generate normalized cluster tensor
    % ------------------------------------------------------------
    DR_j = jj;

    % generate ADPS for cluster
    h = d2d_compute_cluster_tensor(cl, i, grids, quant, cfg, DR_j);
    % cutoff level in relative power
    cutoff_relpwr = 10^(-DR_j/10);

    % ------------------------------------------------------------
    % Step 2: Apply global threshold
    % ------------------------------------------------------------
    % Only retain entries above cutoff
    mask = h >= cutoff_relpwr;

    % If no element survives threshold → skip this DR
    if ~any(mask(:))
        continue;
    end

    % ------------------------------------------------------------
    % Step 3: Scale tensor based on DR
    % ------------------------------------------------------------
    % The scaling enforces that:
    %   max(h2) = cutoff / 10^(-DR/10) * max(h)
    %
    % Since h is normalized (max = 1), scaling factor becomes:
    scale = cutoff / cutoff_relpwr;

    % Apply scaling
    h2 = h * scale;

    % ------------------------------------------------------------
    % Step 4: Compute power mismatch
    % ------------------------------------------------------------
    % Only consider power above threshold
    total_power = sum(h2(mask));

    % Compare with target cluster power
    diff = abs(total_power - target);

    % ------------------------------------------------------------
    % Step 5: Keep best match
    % ------------------------------------------------------------
    if diff < best
        best = diff;
        Hc = h2;
    end

end

end