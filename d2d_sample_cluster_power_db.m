function pwr = d2d_sample_cluster_power_db(delay, d, cfg)
%
% This function computes cluster powers (in dB) based on delay-dependent
% decay models relative to free-space path gain. The power model captures:
%   - exponential decay with delay (log-domain linear)
%   - random fading variation
%   - different behavior for LoS vs NLoS

% ========================================================================
% INPUT
% ========================================================================
% delay : cluster delays [N × 1]
% d     : Tx–Rx distance
% cfg   : configuration (LoS_flag)

% ========================================================================
% OUTPUT
% ========================================================================
% pwr   : cluster powers in dB [N × 1]
% ========================================================================


N = numel(delay);
pwr = zeros(N,1);

for i = 1:N

    % Relative delay in dB domain
    ratio = max(10*log10(delay(i)/d),0);

    if cfg.LoS_flag

        if i == 1
            % LoS component: strong deterministic component
            a1 = 0;    b1 = 7.63;
            a2 = 0;    b2 = 1.85;
        else
            % NLoS clusters in LoS scenario
            a1 = -1.19; b1 = -20.69;
            a2 = -0.17; b2 = 9.26;
        end

    else
        % Pure NLoS scenario
        a1 = -1.81; b1 = -33.41;
        a2 = -0.035; b2 = 10.39;
    end

    % Linear decay + Gaussian(dB-scale) randomness
    pwr(i) = a1*ratio + b1 + (a2*ratio + b2)*randn;

end

end