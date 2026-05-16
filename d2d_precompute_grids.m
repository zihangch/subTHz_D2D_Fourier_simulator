function grids = d2d_precompute_grids(cfg)
% 
% This function precomputes the discrete grids used for the Fourier-domain
% channel representation. The grids define the sampling points for:
%   - delay domain
%   - direction-of-arrival (DoA)
%   - direction-of-departure (DoD)
%
% Separate angular grids are defined for DoA and DoD to support:
%   - asymmetric angular resolutions
%   - different transmitter/receiver array configurations
%   - future extensions (e.g., beamforming, EADF)
%
% All grids are represented as column vectors:
%   - delay: [N_delay × 1]
%   - DoA:   [N_ang_DoA × 1]
%   - DoD:   [N_ang_DoD × 1]
%
% This representation is:
%   - mathematically consistent
%   - easy to manipulate using linear algebra
%   - independent of tensor reshaping logic
%
% Any reshaping required for 3D tensor construction is handled locally
% in downstream functions (e.g., via outer products).

% ========================================================================
% INPUT
% ========================================================================
% cfg : (struct)
%       Configuration structure containing:
%         - N_delay    : number of delay bins
%         - N_ang_DoA  : number of DoA angular bins
%         - N_ang_DoD  : number of DoD angular bins

% ========================================================================
% OUTPUT
% ========================================================================
% grids : (struct)
%         grids.delay : [N_delay × 1] column vector
%                       Discrete delay indices
%
%         grids.DoA   : [N_ang_DoA × 1] column vector
%                       Discrete DoA angles (radians)
%
%         grids.DoD   : [N_ang_DoD × 1] column vector
%                       Discrete DoD angles (radians)
%
%         Note:
%         - Angular grids are uniformly sampled over [0, 2π)
%         - Physical delay = index × delay_res (applied elsewhere)

% ========================================================================

% ------------------------------------------------------------------------
% INPUT VALIDATION
% ------------------------------------------------------------------------
if cfg.N_delay <= 0 || mod(cfg.N_delay,1) ~= 0
    error('cfg.N_delay must be a positive integer.');
end

if cfg.N_ang_DoA <= 0 || mod(cfg.N_ang_DoA,1) ~= 0
    error('cfg.N_ang_DoA must be a positive integer.');
end

if cfg.N_ang_DoD <= 0 || mod(cfg.N_ang_DoD,1) ~= 0
    error('cfg.N_ang_DoD must be a positive integer.');
end

% ------------------------------------------------------------------------
% DELAY GRID
% ------------------------------------------------------------------------
% Discrete delay indices:
%   - Range: 0 to N_delay - 1
%   - Stored as column vector for consistent linear algebra operations
grids.delay = (0:cfg.N_delay-1).' * delay_res;


% ------------------------------------------------------------------------
% DoA GRID (Receiver side)
% ------------------------------------------------------------------------
% Uniform angular sampling over [0, 2π)
% Step size = 2π / N_ang_DoA
% Stored as column vector
grids.DoA = (0:cfg.N_ang_DoA-1).' * (2*pi / cfg.N_ang_DoA);


% ------------------------------------------------------------------------
% DoD GRID (Transmitter side)
% ------------------------------------------------------------------------
% Uniform angular sampling over [0, 2π)
% Step size = 2π / N_ang_DoD
% Stored as column vector
grids.DoD = (0:cfg.N_ang_DoD-1).' * (2*pi / cfg.N_ang_DoD);


end