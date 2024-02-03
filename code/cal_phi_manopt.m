function phi = cal_phi_manopt(D,g)

N = size(g,1);
% Create the problem structure.
manifold = complexcirclefactory(N);
problem.M = manifold;
% Define the problem cost function and its gradient.
problem.cost = @cost;
    function f = cost(phi)
        f = real(phi'*D*phi) + real(phi'*g);
    end
problem.grad = @(phi) problem.M.egrad2rgrad(phi,egrad(phi));
    function grad = egrad(phi)
        grad = 2*D*phi + g;
    end
% Execute the optimization
options.tolgradnorm = 1e-6;
options.maxiter = 100;
options.minstepsize = 1e-6;
options.verbosity = 0;
% checkgradient(problem);
[phi,~,~] = conjugategradient(problem,[],options);

end