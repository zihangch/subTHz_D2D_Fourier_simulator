function [tauD, bA, bD] = d2d_solve_decay_params( ...
    DS, AS_DoA, AS_DoD, DR, quant, cfg)
%
% This function solves the inverse problem of mapping observed spreads
% (after discretization) to true exponential decay parameters.
%
% The process includes:
%   1. Quantization error compensation
%   2. Solving nonlinear inverse mapping using fsolve
%
% The inversion is conditioned on:
%   - delay spread (DS)
%   - angular spreads (AS_DoA, AS_DoD)
%   - effective dynamic range (DR)

% ========================================================================
% INPUT
% ========================================================================
% DS      : observed delay spread
% AS_DoA  : observed angular spread (arrival)
% AS_DoD  : observed angular spread (departure)
% DR      : effective dynamic range (dB)
% quant   : quantization lookup table
% cfg     : configuration structure
%
% ========================================================================
% OUTPUT
% ========================================================================
% tauD   : delay decay constant
% bA     : DoA angular decay constant
% bD     : DoD angular decay constant
% ========================================================================


% ------------------------------------------------------------------------
% STEP 1: SELECT QUANTIZATION TABLE BASED ON DR
% ------------------------------------------------------------------------
% Clamp DR to valid range of lookup table
DR_idx = min(max(round(DR),1), size(quant,1));

if DR_idx > numel(quant)
    error('DR index exceeds quantization table size.');
end

% quantization lookup table
tab = quant{DR_idx,1};


% ------------------------------------------------------------------------
% STEP 2: FIND CLOSEST DISCRETIZED REPRESENTATION
% ------------------------------------------------------------------------
% Table is assumed to store:
%   S_tau_actual, S_tau_theo
%   (extendable to angular spreads if available)

S_tau_actual = tab.S_tau_actual;
S_tau_theo   = tab.S_tau_theo;

% Find closest match in table
[~, min_idx] = min(abs(S_tau_actual - DS), [], 'all', 'linear');
[x,y,z] = ind2sub(size(S_tau_actual), min_idx);


% ------------------------------------------------------------------------
% STEP 3: COMPENSATE QUANTIZATION ERROR
% ------------------------------------------------------------------------
% Remove discretization bias
DS_corr = DS - (S_tau_actual(x,y,z) - S_tau_theo(x,y,z));

% --- Angular spreads (if tables exist) ---
if isfield(tab,'S_AoA_actual') && isfield(tab,'S_AoA_theo')
    ASA_corr = AS_DoA - (tab.S_AoA_actual(x,y,z) - tab.S_AoA_theo(x,y,z));
else
    ASA_corr = AS_DoA;
end

if isfield(tab,'S_AoD_actual') && isfield(tab,'S_AoD_theo')
    ASD_corr = AS_DoD - (tab.S_AoD_actual(x,y,z) - tab.S_AoD_theo(x,y,z));
else
    ASD_corr = AS_DoD;
end


% ------------------------------------------------------------------------
% STEP 4: SOLVE INVERSE PROBLEM
% ------------------------------------------------------------------------
% Solve:
%   f(x, DR, stats) = 0
% where:
%   x = [DS_tmp, bA, bD]

stats = [DS_corr; ASA_corr; ASD_corr];

% Initial guess
x0 = [DS_corr; 0.1; 0.1];

% Nonlinear solve
% "fsolve" solves the problem of f=0
% note: user can also define the epsilon for the residue
opts = optimoptions('fsolve','Display','off');

try
    x_sol = fsolve(@(x) d2d_spread_forward_model_residual(x, DR, stats), x0, opts);
catch
    x_sol = x0;
end

% ------------------------------------------------------------------------
% STEP 5: OUTPUT PARAMETERS
% ------------------------------------------------------------------------
tauD   = x_sol(1);
bA     = x_sol(2);
bD     = x_sol(3);

end