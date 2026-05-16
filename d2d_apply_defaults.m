function cfg = d2d_apply_defaults(cfg)
% 
% This function initializes and completes the configuration structure used
% by the sub-THz D2D Fourier-domain simulator. It assigns default values to
% missing fields while preserving any user-defined parameters.
%
% The function ensures that all required simulation parameters are defined,
% including:
%   - random seed and number of realizations
%   - carrier frequency and scenario type (LoS / NLoS)
%   - Fourier grid resolution (delay domain and angular domains for DoA/DoD)
%   - model parameters such as dynamic range and decay search limits
%   - external data dependencies (e.g., quantization lookup table)
%
% This function does NOT perform:
%   - random sampling (e.g., distances, channels)
%   - simulation execution
%   - parameter validation beyond existence/default assignment
%
% It is intended to provide a clean, deterministic initialization step
% before running the main simulator.

% ========================================================================
% INPUT
% ========================================================================
% cfg : (struct, optional)
%       Configuration structure with user-defined parameters.
%       Missing fields will be automatically filled with default values.
%
%       If empty or not provided, a fully default configuration is created.
%
%       Users may define any subset of fields; unspecified fields will be
%       assigned defaults.

% ========================================================================
% OUTPUT
% ========================================================================
% cfg : (struct)
%       Completed configuration structure with all required fields defined.
%
%       The output cfg includes:
%         - core simulation parameters (seed, N_links, fc, LoS_flag)
%         - Fourier grid parameters (N_delay, N_ang_DoA, N_ang_DoD, delay_res)
%         - model parameters (dynamic_range, max_DR, dom_search_K)
%         - file paths (quant_file)
%
%       This cfg is ready to be used by the main simulator.

% ========================================================================


% If no input is provided, initialize an empty struct
if nargin == 0
    cfg = struct;
end

% ------------------------------------------------------------------------
% CORE SIMULATION SETTINGS
% ------------------------------------------------------------------------

% Base random seed for reproducibility
cfg.seed          = getf(cfg,'seed',102);

% Number of independent channel realizations (Monte Carlo runs)
cfg.N_links       = getf(cfg,'N_links',1);

% Carrier frequency [Hz] (default: sub-THz band)
cfg.fc            = getf(cfg,'fc',145.5e9);

% Scenario type:
%   1 → Line-of-Sight (LoS)
%   0 → Non-Line-of-Sight (NLoS)
cfg.LoS_flag      = getf(cfg,'LoS_flag',0);

% ------------------------------------------------------------------------
% FOURIER GRID PARAMETERS
% ------------------------------------------------------------------------

% Number of delay bins (size of delay-domain grid)
cfg.N_delay       = getf(cfg,'N_delay',1001);

% Number of angular bins for DoA (receiver side)
cfg.N_ang_DoA     = getf(cfg,'N_ang_DoA',36);

% Number of angular bins for DoD (transmitter side)
cfg.N_ang_DoD     = getf(cfg,'N_ang_DoD',36);

% Delay resolution (quantization step in delay domain)
cfg.delay_res     = getf(cfg,'delay_res',0.3);

% ------------------------------------------------------------------------
% MODEL PARAMETERS
% ------------------------------------------------------------------------

% Dynamic range threshold (dB) used for truncation
% This defines the observable range relative to the dominant cluster peak.
% The model is only validated up to 35 dB dynamic range.
cfg.dynamic_range = getf(cfg,'dynamic_range',35);

% --- Validity check: enforce model limit ---
if cfg.dynamic_range > 35
    error('cfg.dynamic_range cannot exceed 35 dB (model validity limit).');
end

if cfg.dynamic_range <= 0
    error('cfg.dynamic_range must be positive.');
end

% Maximum dynamic range index used during power matching
% This should not exceed the global dynamic range.
% We bind it directly to dynamic_range to ensure consistency.
cfg.max_DR = floor(cfg.dynamic_range);

% Number of strongest clusters used in dominant cluster search
% (only relevant for NLoS scenarios)
cfg.dom_search_K = getf(cfg,'dom_search_K',5);

% ------------------------------------------------------------------------
% EXTERNAL DATA FILES
% ------------------------------------------------------------------------

% Quantization error lookup table file
cfg.quant_file    = getf(cfg,'quant_file','quant_err.mat');

end


% ========================================================================
% HELPER FUNCTION
% ========================================================================
function v = getf(s,f,d)
% Returns field value if it exists; otherwise returns default value
if isfield(s,f)
    v = s.(f);
else
    v = d;
end
end