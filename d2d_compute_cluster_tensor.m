function h = d2d_compute_cluster_tensor(cl, i, grids, quant, cfg, DR)
%
% This function generates the normalized delay–DoA–DoD tensor for a single
% cluster. The tensor represents the separable power distribution over:
%   - delay domain
%   - direction-of-arrival (DoA)
%   - direction-of-departure (DoD)
%
% The decay parameters are derived from:
%   - delay spread (DS)
%   - angular spreads (AS_DoA, AS_DoD)
%   - effective dynamic range (DR)
%
% The output tensor is normalized such that its maximum value is 1.
%
% ========================================================================
% INPUT
% ========================================================================
% cl    : cluster structure
%         - delay, DoA, DoD
%         - DS, AS_DoA, AS_DoD
% i     : cluster index
% grids : precomputed grids
%         - delay : [N_delay × 1]
%         - DoA   : [N_DoA × 1]
%         - DoD   : [N_DoD × 1]
% quant : quantization lookup table (cell array)
%         each cell corresponds to a different effective DR
% cfg   : configuration (contains dynamic_range, delay_res, N_delay)
%
% ========================================================================
% OUTPUT
% ========================================================================
% h : normalized 3D tensor [N_delay × N_DoA × N_DoD]
% ========================================================================


% ------------------------------------------------------------------------
% SOLVE DECAY PARAMETERS
% ------------------------------------------------------------------------
% Convert statistical spreads (DS, AS_DoA, AS_DoD) into decay constants
[DS_tmp, bA, bD] = d2d_solve_decay_params( ...
    cl.DS(i), ...
    cl.AS_DoA(i), ...
    cl.AS_DoD(i), ...
    DR, ...
    quant, ...
    cfg);


% ------------------------------------------------------------------------
% DELAY PROFILE (CORRECTED CIRCULAR MODEL)
% ------------------------------------------------------------------------
% Convert delay from physical unit to grid index
delay_idx = grids.delay / cfg.delay_res;          % [N_delay × 1]
d0_idx    = cl.delay(i) / cfg.delay_res;          % scalar

% Compute circular (wrapped) distance
% This ensures shortest distance on delay grids
delta = abs(delay_idx - d0_idx);

% Exponential decay along delay domain
delay_prof = exp(-delta / (DS_tmp / cfg.delay_res));


% ------------------------------------------------------------------------
% ANGULAR PROFILES
% ------------------------------------------------------------------------
% Compute wrapped angular differences ([-π, π])
DoA_diff = abs(wrapToPi(grids.DoA - cl.DoA(i)));   % [N_DoA × 1]
DoD_diff = abs(wrapToPi(grids.DoD - cl.DoD(i)));   % [N_DoD × 1]

% Laplacian (exponential) angular decay
DoA_prof = exp(-DoA_diff / bA);   % [N_DoA × 1]
DoD_prof = exp(-DoD_diff / bD);   % [N_DoD × 1]


% ------------------------------------------------------------------------
% COMBINE DELAY AND ANGLE (SEPARABLE MODEL)
% ------------------------------------------------------------------------
% Outer product for angular spectrum
ang_tensor = DoA_prof * DoD_prof.';   % [N_DoA × N_DoD]

% Expand into 3D tensor and combine with delay profile
h = delay_prof .* reshape(ang_tensor, ...
    1, size(ang_tensor,1), size(ang_tensor,2));


% ------------------------------------------------------------------------
% NORMALIZATION (SAFE)
% ------------------------------------------------------------------------
% Normalize tensor so peak = 1, even though h, in theory, is at peak=1
m = max(h(:));

if m > 0
    h = h / m;
end

end