function W = get_initial_W_SI(C,D,M,N)

% Create the problem structure.
manifold = spherecomplexfactory(M,N);
problem.M = manifold;
% Define the problem cost function and its gradient.

problem.cost = @cost;
    function f = cost(W)
        f = -norm(C*W,'fro')^2 + norm(D*W,'fro')^2;
    end
problem.grad = @(W) problem.M.egrad2rgrad(W,egrad(W));
    function g = egrad(W)
        g = -2*(C'*C)*W + 2*(D'*D)*W;
    end

% Execute the optimization
options.tolgradnorm = 1e-6;
options.maxiter = 1000;
options.minstepsize = 1e-6;
options.verbosity = 0;
% checkgradient(problem);

[W,~,~] = conjugategradient(problem,[],options);


end