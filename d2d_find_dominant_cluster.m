function idx = d2d_find_dominant_cluster(cl, grids, quant, cfg)
% 
% This function identifies the dominant cluster in a given channel
% realization. The dominant cluster is defined as the one that produces
% the largest observable peak in the joint delay–DoA–DoD spectrum under
% a fixed dynamic range (DR) constraint.
%
% Algorithm:
%   - LoS scenario:
%       → The first cluster (direct path) is always dominant
%
%   - NLoS scenario:
%       → Select top-K strongest clusters (by large-scale power)
%       → For each candidate:
%           1. Construct normalized cluster tensor
%           2. Apply DR-based thresholding
%           3. Scale tensor to match cluster power
%           4. Evaluate resulting peak
%       → Choose cluster with maximum peak
%
% This ensures a fair comparison across clusters by accounting for:
%   - decay behavior (delay + angular spreads)
%   - discretization effects
%   - power normalization
%
% ========================================================================
% INPUT
% ========================================================================
% cl    : cluster structure
%         - pwr      : cluster powers (dB)
%         - DS       : delay spreads
%         - AS_DoA   : angular spread (DoA)
%         - AS_DoD   : angular spread (DoD)
%         - delay    : cluster delays
%         - DoA/DoD  : angular centroids
%
% grids : precomputed grids (delay, DoA, DoD)
%
% quant : quantization lookup table (cell array)
%         each cell corresponds to a different effective dynamic range
%
% cfg   : configuration structure
%         - dynamic_range : global DR (dB)
%         - dom_search_K  : number of candidate clusters
%
% ========================================================================
% OUTPUT
% ========================================================================
% idx : index of dominant cluster
% ========================================================================


% ------------------------------------------------------------------------
% INPUT VALIDATION
% ------------------------------------------------------------------------
% Ensure quantization table is correctly loaded
if ~iscell(quant)
    error('quantization lookup table must be provided as a cell array.');
end


% ------------------------------------------------------------------------
% LoS SCENARIO: direct path always dominant
% ------------------------------------------------------------------------
if cfg.LoS_flag
    idx = 1;
    return;
end


% ------------------------------------------------------------------------
% NLoS SCENARIO: select top-K strongest clusters
% ------------------------------------------------------------------------

% Number of candidates (cfg already provides default)
K = min(cfg.dom_search_K, cl.N);

% Select top-K clusters based on large-scale power
[~, cand] = maxk(cl.pwr, K);

% Safety check (should not happen, but prevents crash)
if isempty(cand)
    idx = 1;
    return;
end


% ------------------------------------------------------------------------
% GLOBAL DR THRESHOLD
% ------------------------------------------------------------------------
% Convert DR (dB) to linear scale threshold ratio
DR = cfg.dynamic_range;
cutoff_ratio = 10^(-DR/10);

% Note:
% h (cluster tensor) is assumed normalized such that max(h) = 1.
% Therefore, cutoff_ratio directly represents the observable threshold.


% ------------------------------------------------------------------------
% DOMINANT CLUSTER SEARCH
% ------------------------------------------------------------------------
best_peak = -inf;
idx = cand(1);   % initialize with strongest cluster (by power)

for i = cand

    % ------------------------------------------------------------
    % Step 1: construct normalized cluster tensor
    % ------------------------------------------------------------
    DR = cfg.dynamic_range;

    h = d2d_compute_cluster_tensor(cl, i, grids, quant, cfg, DR);

    % ------------------------------------------------------------
    % Step 2: apply DR threshold (retain only observable region)
    % ------------------------------------------------------------
    mask = h >= cutoff_ratio;

    % Skip cluster if nothing survives threshold
    if ~any(mask(:))
        continue;
    end

    % ------------------------------------------------------------
    % Step 3: scale tensor to match cluster power
    % ------------------------------------------------------------
    % Convert cluster power from dB to linear scale
    P_lin = 10^(cl.pwr(i)/10);

    % Extract surviving region after thresholding
    h_masked = h .* mask;

    % Compute total power (discrete approximation)
    % Note: sum() approximates integral over delay–angle grid
    total_pwr = sum(h_masked(:));

    % Skip degenerate case
    if total_pwr == 0
        continue;
    end

    % Normalize and scale to match cluster power
    h_scaled = h_masked / total_pwr * P_lin;

    % ------------------------------------------------------------
    % Step 4: evaluate resulting peak
    % ------------------------------------------------------------
    peak_val = max(h_scaled(:));

    % Update dominant cluster
    if peak_val > best_peak
        best_peak = peak_val;
        idx = i;
    end

end

end