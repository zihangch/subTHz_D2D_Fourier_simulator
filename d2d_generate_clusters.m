function cl = d2d_generate_clusters(cfg, d)
%
% This function generates stochastic channel clusters for a sub-THz D2D
% channel model. Each cluster is characterized by:
%   - delay (arrival time)
%   - angular centroids (DoA, DoD)
%   - delay spread (DS)
%   - angular spreads (AS_DoA, AS_DoD)
%   - power (large-scale fading)
%
% The model supports both LoS and NLoS scenarios with different statistics.

% ========================================================================
% INPUT
% ========================================================================
% cfg : configuration struct (LoS_flag, delay_res, N_delay, etc.)
% d   : Tx–Rx distance [m]

% ========================================================================
% OUTPUT
% ========================================================================
% cl : struct containing cluster parameters
% ========================================================================


% ------------------------------------------------------------------------
% NUMBER OF CLUSTERS
% ------------------------------------------------------------------------
% number of clusters are rounded to nearest integer
if cfg.LoS_flag
    N = 1 + round(gamrnd(3.39,3.35));
else
    N = 1 + round(gamrnd(1.89,46.72));
end


% ------------------------------------------------------------------------
% CLUSTER ARRIVAL TIMES (EXCESS DELAY AGAINST LOS DISTANCE)
% ------------------------------------------------------------------------
if cfg.LoS_flag
    tau = [0; sort(gamrnd(0.86,32.40,[max(N-1,0),1]))];
else
    tau = sort(gamrnd(1.51,35,[N,1]));
end

% Enforce delay support within FFT grid
max_delay = (cfg.N_delay-1) * cfg.delay_res;
idx = tau > max_delay;

% Resample exceeding delays to guarantee they're within valid supports
while ~isempty(idx)
    if cfg.LoS_flag
        tau(idx) = gamrnd(0.86,32.40,[length(idx),1]);
    else
        tau(idx) = gamrnd(1.51,35,[length(idx),1]);
    end
    tau = sort(tau,'ascend');
    idx = tau > max_delay;
end

% Quantize to nearest delay grid (FFT-compatible discrete taps)
delay = round((d + tau)/cfg.delay_res)*cfg.delay_res;


% ------------------------------------------------------------------------
% ANGULAR CENTROIDS (DoA / DoD)
% ------------------------------------------------------------------------
DoA = zeros(N,1);
DoD = zeros(N,1);

if cfg.LoS_flag

    % Direct LoS cluster
    DoA(1) = 0;
    DoD(1) = 0;

    % Other clusters follow bi-modal Laplacian
    for ii = 2:N

        % -------- DoA --------
        if rand < 0.7344
            dev = exprnd(0.292);
            sign = (rand > 0.5)*2 - 1;
            DoA(ii) = wrapToPi(sign * dev);
        else
            dev = exprnd(0.068);
            sign = (rand > 0.5)*2 - 1;
            DoA(ii) = wrapToPi(pi + sign*dev);
        end

        % -------- DoD --------
        if rand < 0.318
            dev = exprnd(0.151);
            sign = (rand > 0.5)*2 - 1;
            DoD(ii) = wrapToPi(sign * dev);
        else
            dev = exprnd(1.988);
            sign = (rand > 0.5)*2 - 1;
            DoD(ii) = wrapToPi(pi + sign*dev);
        end

    end

else
    % NLoS: isotropic scattering (wrapped to [-π, π))
    DoA = wrapToPi(2*pi*rand(N,1));
    DoD = wrapToPi(2*pi*rand(N,1));
end


% ------------------------------------------------------------------------
% DELAY SPREAD (DS)
% ------------------------------------------------------------------------
if cfg.LoS_flag
    DS = [lognrnd(-1.42,0.43); lognrnd(-1.42,0.52,[max(N-1,0),1])];
else
    DS = lognrnd(-0.79,0.89,[N,1]);
end


% ------------------------------------------------------------------------
% ANGULAR SPREADS
% ------------------------------------------------------------------------
if cfg.LoS_flag
    AS_DoA = [lognrnd(-1.98,0.042); lognrnd(-2.30,0.29,[max(N-1,0),1])];
    AS_DoD = [lognrnd(-1.96,0.048); lognrnd(-2.26,0.32,[max(N-1,0),1])];
else
    AS_DoA = lognrnd(-2.18,0.35,[N,1]);
    AS_DoD = lognrnd(-2.17,0.38,[N,1]);
end


% ------------------------------------------------------------------------
% CLUSTER POWER
% ------------------------------------------------------------------------
pwr = d2d_sample_cluster_power_db(delay, d, cfg);


% ------------------------------------------------------------------------
% OUTPUT STRUCT
% ------------------------------------------------------------------------
cl = struct('N',N,'delay',delay,'DoA',DoA,'DoD',DoD,...
            'DS',DS,'AS_DoA',AS_DoA,'AS_DoD',AS_DoD,...
            'pwr',pwr);

end