% Obtain an initial phi by maximizing the channel gains.
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Download this paper at: https://ieeexplore.ieee.org/document/10364735
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-28
% Inputs: A, b
% Outputs: phi: RIS reflecting coefficients

function phi = get_initial_phi(A,b)

[N,~] = size(b);

% Create the problem structure.
manifold = complexcirclefactory(N);
problem.M = manifold;
% Define the problem cost function and its gradient.

problem.cost = @cost;
    function f = cost(x)
        f = -x'*A*x - real(x'*b);
    end
problem.grad = @(x) problem.M.egrad2rgrad(x,egrad(x));
    function g = egrad(x)
        g = -2*A*x - b;
    end

% Execute the optimization
options.tolgradnorm = 1e-6;
options.maxiter = 1000;
options.minstepsize = 1e-6;
options.verbosity = 0;
% checkgradient(problem);
[phi,~,~] = conjugategradient(problem,[],options);

end


