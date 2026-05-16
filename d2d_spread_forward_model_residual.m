function F = d2d_spread_forward_model_residual(x, DR, stats)
%
% This function evaluates the discrepancy between:
%   - the observed (target) delay and angular spreads
%   - the theoretically predicted spreads derived from decay parameters
%
% The underlying model describes the analytical relationship between:
%   - exponential decay in delay domain
%   - Laplacian decay in angular domains (DoA / DoD)
%   - truncation induced by a finite dynamic range (DR)
%
% Given decay parameters and DR, the function computes:
%   - delay spread
%   - DoA angular spread
%   - DoD angular spread
%
% These predicted spreads are compared against the desired statistics,
% and the residual is returned for use in a nonlinear solver (e.g., fsolve).
%
% In summary:
%   Forward model:
%       (decay constants, DR) → (observed spreads)
%
%   Inverse problem (via solver):
%       (observed spreads, DR) → (decay constants)
%
% ========================================================================
% INPUT
% ========================================================================
% x     : [3×1] decay parameters
%         x(1) → delay exponential decay constant
%         x(2) → DoA angular decay constant
%         x(3) → DoD angular decay constant
%
% DR    : effective dynamic range (dB)
%
% stats : [3×1] target spreads
%         stats(1) → delay spread
%         stats(2) → DoA angular spread
%         stats(3) → DoD angular spread
%
% ========================================================================
% OUTPUT
% ========================================================================
% F : [3×1] residual vector
%     predicted_spreads(x, DR) - stats
% ========================================================================


% ------------------------------------------------------------------------
% PARAMETER REORDERING
% ------------------------------------------------------------------------
% Ensure b ≤ c for consistent piecewise formulation
a = x(1);                  % delay decay constant
b = min(x(2), x(3));       % smaller angular decay constant
c = max(x(2), x(3));       % larger angular decay constant


% ------------------------------------------------------------------------
% DYNAMIC RANGE NORMALIZATION
% ------------------------------------------------------------------------
% Convert DR (dB) to exponential truncation parameter
t = DR / 10 * log(10);


% ------------------------------------------------------------------------
% REGION DEFINITIONS
% ------------------------------------------------------------------------
% Define truncation boundaries in delay-angle domain
W2 = a * (t - pi/b - pi/c);
Wb = a * (t - pi/b);
Wc = a * (t - pi/c);
W  = a * t;


% ------------------------------------------------------------------------
% MOMENT COMPUTATION UNDER TRUNCATION
% ------------------------------------------------------------------------
% These quantities are derived analytically from truncated integrals:
%   P_bar : total power after truncation
%   Ex    : first moment (mean delay)
%   Ex2   : second moment (delay variance)
%   Ey    : angular moment (DoA)
%   Ez    : angular moment (DoD)

if b <= c

    % ------------------------------------------------------------
    % Case 1: No angular truncation (DR dominates)
    % ------------------------------------------------------------
    if Wc <= 0

        P_bar = 1 - exp(-t) - exp(-t)*(t + t^2/2);

        Ex = a * (1 - exp(-t)*(1+t) - exp(-t)*(t^2/2 + t^3/6));

        Ex2 = a^2 * (2 - exp(-t)*(t^2 + 2*t + 2) ...
                       - exp(-t)*(t^3/3 + t^4/12));

        Ey = 1/(b^2+1) - exp(-t)/b^2 ...
            + exp(-t)/(b^2*(b^2+1))*(cos(b*t)-b*sin(b*t));

        Ez = 1/(c^2+1) - exp(-t)/c^2 ...
            + exp(-t)/(c^2*(c^2+1))*(cos(c*t)-c*sin(c*t));


    % ------------------------------------------------------------
    % Case 2: Partial angular truncation (one dimension)
    % ------------------------------------------------------------
    elseif Wc > 0 && Wb <= 0

        P_bar = 1 - exp(-pi/c) - exp(-t)*pi/c*(1+t) ...
                + exp(-t)/2*pi^2/c^2;

        Ex = a * (1 - exp(-pi/c) ...
                - exp(-t)/2*pi/c*(t^2+2*t+2) ...
                + exp(-t)/2*pi^2/c^2*(t+1) ...
                - exp(-t)/6*pi^3/c^3);

        Ex2 = a^2 * (2*(1-exp(-pi/c)) ...
                - exp(-t)*pi/c*(2+2*t+t^2+t^3/3) ...
                + exp(-t)*pi^2/c^2*(1+t+t^2/2) ...
                - exp(-t)/3*pi^3/c^3*(1+t) ...
                + exp(-t)/12*pi^4/c^4);

        Ey = (1-exp(-pi/c))/(b^2+1) ...
            + exp(-t)/(b^2*(b^2+1)) ...
              *(cos(b*t)-cos(b*(t-pi/c)) ...
              - b*sin(b*t)+b*sin(b*(t-pi/c)));

        Ez = (1+exp(-pi/c))/(c^2+1) - 2*exp(-t)/c^2;


    % ------------------------------------------------------------
    % Case 3: Both angular domains partially truncated
    % ------------------------------------------------------------
    elseif Wb > 0 && W2 <= 0

        P_bar = 1 - exp(-pi/b) - exp(-pi/c) ...
              + exp(-t)*(1+t)*(1-pi/b-pi/c) ...
              + exp(-t)/2*(t^2+pi^2/b^2+pi^2/c^2);

        Ex = a * (1 - exp(-pi/b) - exp(-pi/c) ...
            + exp(-t)*(1+t+t^2/2)*(1-pi/b-pi/c) ...
            + exp(-t)/2*(1+t)*(pi^2/b^2+pi^2/c^2) ...
            + exp(-t)/6*(t^3-pi^3/b^3-pi^3/c^3));

        Ex2 = a^2 * (2*(1-exp(-pi/b)-exp(-pi/c)) ...
            + exp(-t)*(1+t+t^2/2)*(2-2*pi/b-2*pi/c ...
            + pi^2/b^2+pi^2/c^2) ...
            + exp(-t)/3*(t^3*(1-pi/b-pi/c) ...
            - (1+t)*(pi^3/b^3+pi^3/c^3)) ...
            + exp(-t)/12*(t^4+pi^4/b^4+pi^4/c^4));

        Ey = (1+exp(-pi/b)-exp(-pi/c))/(b^2+1) ...
            - exp(-t)/b^2 ...
            - exp(-t)/(b^2*(b^2+1)) ...
              *(cos(b*pi/c-b*t)+b*sin(b*pi/c-b*t));

        Ez = (1+exp(-pi/c)-exp(-pi/b))/(c^2+1) ...
            - exp(-t)/c^2 ...
            - exp(-t)/(c^2*(c^2+1)) ...
              *(cos(c*pi/b-c*t)+c*sin(c*pi/b-c*t));


    % ------------------------------------------------------------
    % Case 4: Strong truncation (fully bounded region)
    % ------------------------------------------------------------
    elseif W2 > 0

        P_bar = (1-exp(-pi/b))*(1-exp(-pi/c)) ...
              - pi^2/(b*c)*exp(-t);

        Ex = a * ((1-exp(-pi/b))*(1-exp(-pi/c)) ...
            + exp(-t)/(2*b*c)*(pi^3/b+pi^3/c-2*(1+t)*pi^2));

        Ex2 = a^2 * (2*(1-exp(-pi/b))*(1-exp(-pi/c)) ...
            - exp(-t)/(b*c)*pi^2*(t^2+2*t+2) ...
            + exp(-t)/(b*c)*pi^2*(pi/b+pi/c)*(1+t) ...
            - exp(-t)/(b*c)*pi^2*(pi^2/(3*b^2)+pi^2/(2*b*c)+pi^2/(3*c^2)));

        Ey = (1+exp(-pi/b))*(1-exp(-pi/c))/(b^2+1);
        Ez = (1+exp(-pi/c))*(1-exp(-pi/b))/(c^2+1);

    end
end


% ------------------------------------------------------------------------
% FINAL SPREAD COMPUTATION
% ------------------------------------------------------------------------
% Convert moments to RMS spreads and compute residual

if stats(2) <= stats(3)
    F = [sqrt(Ex2/P_bar - (Ex/P_bar)^2); ...
         sqrt(1 - (Ey/P_bar)^2); ...
         sqrt(1 - (Ez/P_bar)^2)] - stats;
else
    % Swap angular components if needed
    F = [sqrt(Ex2/P_bar - (Ex/P_bar)^2); ...
         sqrt(1 - (Ez/P_bar)^2); ...
         sqrt(1 - (Ey/P_bar)^2)] - stats;
end

end