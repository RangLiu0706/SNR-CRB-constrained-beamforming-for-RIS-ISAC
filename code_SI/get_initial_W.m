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

[W,aa,bb] = conjugategradient(problem,[],options);


end