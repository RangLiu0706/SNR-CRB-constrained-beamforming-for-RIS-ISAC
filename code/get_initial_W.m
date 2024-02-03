% Obtain an initial W by maximizing the received signal power.
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: C, M: number of antennas, N: number of RIS elements
% Outputs: W: transmit beamforming

function W = get_initial_W(C,M,N)

% Create the problem structure.
manifold = spherecomplexfactory(M,N);
problem.M = manifold;
% Define the problem cost function and its gradient.

problem.cost = @cost;
    function f = cost(W)
        f = -norm(C*W,'fro')^2;
    end
problem.grad = @(W) problem.M.egrad2rgrad(W,egrad(W));
    function g = egrad(W)
        g = -2*(C'*C)*W;
    end

% Execute the optimization
options.tolgradnorm = 1e-6;
options.maxiter = 1000;
options.minstepsize = 1e-6;
options.verbosity = 0;
% checkgradient(problem);

[W,~,~] = conjugategradient(problem,[],options);


end