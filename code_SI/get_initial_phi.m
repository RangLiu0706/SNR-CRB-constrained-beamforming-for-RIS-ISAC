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
[phi,aa,bb] = conjugategradient(problem,[],options);

end


