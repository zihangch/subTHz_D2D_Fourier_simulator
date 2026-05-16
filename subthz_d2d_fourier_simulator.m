function sim = subthz_d2d_fourier_simulator(cfg)
% 
% This function generates a sub-THz D2D channel using a Fourier-domain 
% cluster-based model. Each channel realization is constructed as a 3D 
% tensor over:
%   - delay domain
%   - direction-of-arrival (DoA)
%   - direction-of-departure (DoD)
%
% The simulator includes:
%   - stochastic cluster generation (delays, angles, spreads, powers)
%   - quantization-aware spread correction (via lookup table)
%   - dominant cluster detection (to define dynamic range)
%   - threshold-based truncation (emulating measurement/DR limits)
%   - power-consistent tensor reconstruction
%
% Each realization is independent and reproducible via deterministic seeding.
% The model is suitable for high-resolution channel modeling and analysis
% in sub-THz systems.

% ========================================================================
% INPUT
% ========================================================================
% cfg : (struct, optional)
%       Configuration structure controlling simulation behavior.
%       If empty or not provided, default parameters are applied.
%
%       Key fields (all optional unless stated):
%         - seed              : base random seed (default: 102)
%         - N_links           : number of channel realizations (default: 1)
%         - fc                : carrier frequency [Hz]
%         - LoS_flag          : 1 = LoS, 0 = NLoS
%         - distance          : scalar Tx–Rx distance [m] (optional)
%         - N_delay           : number of delay bins
%         - N_ang_DoA         : number of DoA angular bins
%         - N_ang_DoD         : number of DoD angular bins
%         - delay_res         : delay resolution
%         - dynamic_range     : threshold in dB (e.g., 35)
%         - quant_file        : path to quantization error table
%
%       Note:
%       - distance must be scalar if provided
%       - multiple realizations are handled via N_links

% ========================================================================
% OUTPUT
% ========================================================================
% sim : struct array of length N_links
%       Each element contains:
%
%       sim(nn).H
%           3D channel tensor of size:
%           [N_delay × N_ang_DoA × N_ang_DoD]
%           representing delay–DoA–DoD response
%
%       sim(nn).meta
%           struct with metadata:
%             - distance           : Tx–Rx distance [m]
%             - N_cluster          : number of clusters
%             - dominant_cluster   : index of dominant cluster
%             - path_gain          : Friis path gain (linear)
%
% ========================================================================

% ---------------- Input handling ----------------
if nargin == 0 || isempty(cfg)
    % empty input, configuration takes all default values
    cfg = d2d_apply_defaults();
else
    % all non-assigned fields are assigned with default values
    cfg = d2d_apply_defaults(cfg);
end

% ---------------- Initialization ----------------
% random seed initialization
rng(cfg.seed);  
% load quantization lookup table
tmp = load(cfg.quant_file);
quant = tmp.quant_err;   % extract cell array
% compute Fourier grids
grids = d2d_precompute_grids(cfg);

% Total number of channel realizations to generate (each produces one H tensor)
N = cfg.N_links;
sim = repmat(struct('H',[],'meta',[]), N, 1);

% ================================================================
% DISTANCE HANDLING (INITIALIZATION STAGE)
% ================================================================
% Supported cases:
%   1) cfg.distance is scalar
%        → same distance for all realizations
%   2) cfg.distance is vector with length N_links
%        → one distance per realization
%   3) cfg.distance not provided
%        → randomly generate distances in [1, 100] m

if isfield(cfg,'distance')

    if numel(cfg.distance) == 1
        % --- Single distance for all realizations ---
        d_vec = repmat(cfg.distance, N, 1);

    elseif numel(cfg.distance) == N
        % --- User-provided distance per realization ---
        d_vec = cfg.distance(:);

    else
        error(['cfg.distance must be either a scalar or a vector of ' ...
               'length cfg.N_links.']);
    end

    % --- Validity check: model range ---
    if any(d_vec < 1 | d_vec > 100)
        error('All distance values must be within [1, 100] meters.');
    end

else
    % --- Default: random distances in valid range ---
    d_vec = 1 + (100-1)*rand(N,1);
end

% ---------------- Main loop ----------------
for nn = 1:N

    % ============================================================
    % RANDOM SEED HANDLING
    % ============================================================
    % Two supported modes:
    %   1) cfg.seed is scalar:
    %        → use (seed + nn) to generate independent realizations
    %   2) cfg.seed is a vector of length N_links:
    %        → use user-provided seed for each realization
    
    if numel(cfg.seed) == 1
        % --- Scalar seed: deterministic progression ---
        rng(cfg.seed + nn);
        
    elseif numel(cfg.seed) == N
        % --- Vector seeds: one seed per realization ---
        rng(cfg.seed(nn));
        
    else
        error(['cfg.seed must be either a scalar or a vector of length ' ...
            'cfg.N_links.']);
    end

    % --- Distance (precomputed) ---
    d = d_vec(nn);

    % --- Friis path gain ---
    PG = (3e8/(4*pi*d*cfg.fc))^2;

    % --- Generate cluster parameters ---
    cl = d2d_generate_clusters(cfg, d);

    % --- Determine dominant cluster ---
    dom_idx = d2d_find_dominant_cluster(cl, grids, quant, cfg);

    % --- Threshold ---
    DR = cfg.dynamic_range;
    h_dom = d2d_compute_cluster_tensor(cl, dom_idx, grids, quant, cfg, DR);
    % cutoff level, in relative power
    cutoff = max(h_dom(:)) * 10^(-cfg.dynamic_range/10);
    % take out all surviving grids
    mask = h_dom >= cutoff;
    h_masked = h_dom .* mask;
    % dominant cluster peak power = 1
    % normalized dominant total cluster power, in absolute power
    h_scale = h_masked / sum(h_masked(:)) * 10^((cl.pwr(dom_idx)+10*log10(PG))/10);
    % cutoff level, in absolute power
    cutoff_global = max(h_scale(:)) * 10^(-cfg.dynamic_range/10);

    % --- Channel ---
    H = zeros(cfg.N_delay, cfg.N_ang_DoA, cfg.N_ang_DoD);

    for ii = 1:cl.N
        Hc = d2d_optimize_cluster_power(cl, ii, grids, quant, cutoff_global, PG, cfg);
        if ~isempty(Hc)
            H = H + Hc;
        end
    end

    % --- Store ---
    sim(nn).H = H;
    sim(nn).meta = struct( ...
        'distance', d, ...
        'N_cluster', cl.N, ...
        'dominant_cluster', dom_idx, ...
        'path_gain', PG );

end

end