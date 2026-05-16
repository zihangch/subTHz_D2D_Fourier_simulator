1. Overview

This repository contains a stochastic Fourier-domain channel simulator for sub-THz outdoor D2D communications.
It models multipath propagation using a cluster-based representation in the delay–angle domain, where each cluster
is described by physically meaningful delay and angular spreads.

The simulator generates high-resolution delay–DoA–DoD channel tensors, incorporating:
  - stochastic cluster generation (delay, angle, spreads, power)
  - quantization-aware spread correction (via lookup table)
  - dynamic-range (DR) dependent truncation effects
  - inversion-based mapping from observed spreads to decay parameters
  - consistent power normalization under DR constraints

The resulting channel tensors are suitable for high-resolution channel analysis, beamforming studies, and
model validation in sub-THz systems.


2. What is needed to run

- MATLAB (Optimization Toolbox required for fsolve)
- All .m files in this folder added to the MATLAB path
- quant_err.mat file (quantization error lookup table)


3. Quick Example

cfg = struct('distance', 30, 'LoS_flag', 1);
sim = subthz_d2d_fourier_simulator(cfg);

imagesc(abs(sim(1).H(:,:,1)))
title('Channel magnitude (delay × DoA)')


4. Structure

/simulator
│
├── subthz_d2d_fourier_simulator.m   % Main simulation script
│
├── Core functions
│   ├── d2d_generate_clusters.m        % Generate cluster parameters
│   ├── d2d_compute_cluster_tensor.m   % Build DR-dependent cluster tensor
│   ├── d2d_find_dominant_cluster.m    % Identify dominant cluster
│   ├── d2d_optimize_cluster_power.m   % Match cluster power under DR constraint
│
├── Inversion / modeling
│   ├── d2d_solve_decay_params.m       % Solve decay parameters via inversion
│   ├── d2d_spread_forward_model_residual.m
│       % Forward model: decay parameters → observed spreads
│
├── Statistical sampling
│   ├── d2d_sample_cluster_power_db.m  % Delay-dependent cluster power model
│
├── Geometry & grids
│   ├── d2d_precompute_grids.m         % Delay and angular grids
│
├── Configuration
│   ├── d2d_apply_defaults.m           % Default configuration setup
│
├── Data
│   ├── quant_err.mat                  % Quantization error lookup table


5. Key Modeling Features

- Delay–angle separable channel representation
- DR-dependent tensor truncation (critical for realistic observation)
- Inversion-based parameter estimation (spreads → decay constants)
- Dominant cluster selection based on observable peak (not total power)
- Consistent treatment of discretization and quantization effects


6. Notes

- The simulator operates on discrete grids; delay and angular resolutions
  must be chosen carefully for accuracy.
- The quantization lookup table (quant_err.mat) must match the assumed grid setup.
- Dynamic range (DR) plays a central role and affects both tensor shape and power.


7. Output Format

sim(nn).H:
    3D channel tensor of size:
        [N_delay × N_ang_DoA × N_ang_DoD]

sim(nn).meta:
    distance            : Tx–Rx distance (m)
    N_cluster           : number of clusters
    dominant_cluster    : index of dominant cluster
    path_gain           : Friis path gain (linear)

8. Quantization Lookup Table (quant_err.mat)

The file `quant_err.mat` contains a precomputed lookup table used for
quantization-aware spread correction in the simulator. It compensates for
discretization effects introduced by the finite delay–angle grid.

Structure:
----------
- The variable `quant_err` is a cell array of size [N_DR × 1], where:
    N_DR ≈ 35 (corresponding to dynamic range values from 1 dB to 35 dB)

- Each cell corresponds to a specific effective dynamic range (DR) and contains
  a struct with the following fields:

    S_tau_actual   : [50 × 23 × 23] double
    S_phiR_actual  : [50 × 23 × 23] double
    S_phiT_actual  : [50 × 23 × 23] double

    S_tau_theo     : [50 × 23 × 23] double
    S_phiR_theo    : [50 × 23 × 23] double
    S_phiT_theo    : [50 × 23 × 23] double

Meaning:
--------
- "actual" refers to observed spreads after discretization on the delay–angle grid
- "theo" refers to the corresponding theoretical spreads (continuous-domain)

The lookup table enables inversion from observed spreads to underlying
decay parameters by compensating discretization bias.

Grid Definition:
----------------
The 3D matrices are indexed by:

1) Delay decay constant (τ-domain):
   - Size: 50
   - Range: 0.1 : 0.1 : 5.0

2) Angular decay constant (DoA):
   - Size: 23
   - Range: 0.05 : 0.02 : 0.49

3) Angular decay constant (DoD):
   - Size: 23
   - Same range as DoA

Thus, each entry in the table corresponds to a triplet:
    (delay decay, DoA decay, DoD decay)

Usage in Simulator:
-------------------
- Given observed spreads (DS, AS_DoA, AS_DoD), the simulator:
    1. Finds the closest discretized entry in S_*_actual
    2. Computes bias relative to S_*_theo
    3. Compensates the observed spreads
    4. Solves for true decay parameters via nonlinear inversion

Important Notes:
----------------
- The lookup table is DR-dependent: different cells correspond to different
  truncation levels of the channel (i.e., observable dynamic range).
- The table must be consistent with:
    - delay resolution (cfg.delay_res)
    - angular grid size (cfg.N_ang_DoA, cfg.N_ang_DoD)
- Mismatch between lookup table and simulator grid may lead to incorrect
  spread inversion.
