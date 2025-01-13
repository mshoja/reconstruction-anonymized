function [costfun] = Create_Hermite_Cost_Parametric(data, dt, J, K, mu, sigma, npars, refined, analytic_diffs, fraction)
    % we create an efficient costfunction and the used mu and sd functions
    mesh_fine = linspace(min(data), max(data), 2000);
    if analytic_diffs
        [mufuns, nder_mu] = funderivatives(mu, npars, 2 * K - 2);
        [sigmafuns, nder_sigma] = funderivatives(sigma, npars, 2 * K - 2);
    else
        mesh_fine1 = linspace(min(data), max(data), 40000);
        nder_mu = 2 * K - 2;
        nder_sigma = 2 * K - 2;
        mufuns = @(x,par)numeric_diff(mu,x, par, 2 * K - 2, mesh_fine1);
        sigmafuns = @(x,par)numeric_diff(sigma,x,par, 2 * K - 2,mesh_fine1);
    end
 
    
    %if we have a different number of derivatives for mu and sigma and we
    %have not additive noise (=nder_sigma>0) then we use the maximum of
    %both. This saves creating a lot of different EZ functions if we
    %combine different kinds of splines for mu and sigma. 
    if nder_mu ~= nder_sigma && nder_sigma > 0
        nder_mu = max(nder_mu, nder_sigma);
        nder_sigma = max(nder_mu, nder_sigma);
    end
    if nder_mu < 2 * K - 2
        mufuns = funderivatives(mu, npars, nder_mu);
    end
    if nder_sigma < 2 * K - 2
        sigmafuns = funderivatives(sigma, npars, nder_sigma);
    end
    EZ = shared.get_hermite(J, K, nder_mu, nder_sigma, refined);

    if isempty(fraction)
        N = length(data);
        X0 = data(1:end - 1);
        X = data(2:end);
    else
        sample_pairs = readmatrix('sample_pairs.csv');
        X0 = sample_pairs(1:2:end);
        X = sample_pairs(2:2:end);
        N=length(X0)+1;
    end
    if size(X0, 1) > 1
        X0 = X0.';
        X = X.';
    end
    [H, H_coeff] = shared.Hfun(J,true);
    %mesh_fine is used to approximate the integral of 1/sigma
    %    mesh_fine = linspace(min(data), max(data), 2000);
    iinvsigma = intsigma(sigma, npars);
    costfun = @(par)Cost_Parametric(par, J, N, X0,X, dt,EZ, H, H_coeff, shared.Coefffun(J),mesh_fine,iinvsigma ,mufuns,sigmafuns,sigma);
end
% 
% function res = normfun(fun, x, par)
%     cpar = num2cell(par);
%     res = fun(x, cpar{:}) + zeros(size(x));
% end
function intfunc = intsigma(fun, npars)
    pars = sym('par', [1, npars]);
    syms x F x0 dt
    F = fun(x, pars);
    Z = int(1 / F, x);
    Z = (Z - subs(Z, x, x0));
    Z = simplify(Z);
    intfunc = makevectfun(matlabFunction(Z, 'Vars', [x0 x pars]));
end

function res = numeric_diff(fun, x, pars, nder, mesh_fine)
    F = fun(x, pars);
    res = cell(nder + 1, 1);
    res{1} = F;
    dx = mesh_fine(2) - mesh_fine(1);
    Fmesh = fun(mesh_fine, pars);
    for i = 1:nder
        Fmesh = [diff(Fmesh) / dx NaN];
        if all(Fmesh(~isnan(Fmesh))<1E-8)
            %if the derivatives are almost zero we need to set them to
            %zero, else the errors can blow up in higher derivatives
            Fmesh=zeros(size(Fmesh));
        end 
        res{i + 1} = interp1(mesh_fine, Fmesh, x);
    end
end

function [res, nnonzero] = funderivatives(fun, npars, nder)
    pars = sym('par', [1, npars]);
    syms x F
    nnonzero = nder;
    F = fun(x, pars);
    res = cell(nder + 1, 1);
    for i = 1:nder + 1
        if F == 0
            res{i} = 0;
            nnonzero = nnonzero - 1;
        else
            res{i} = makevectfun(matlabFunction(F, 'vars', [x pars]));
            F = diff(F, x);
        end
    end
    
end

function fun = makevectfun(fun)
    f = func2str(fun);
    %remove @(x,par1,par2 etc)
    [a, b] = regexp(f, '(par[0-9]*[,)])*', 'start', 'end', 'once');
    f = [f(1:a - 1) 'par' f(b:end)];
    %get all variable names
    vars = regexp(f, '\<[a-zA-Z][a-zA-Z_0-9]*', 'match');
    %is there no x? add +zeros(size(x)) to the function
    if ~any(strcmp(vars, 'x'))
        f = [f '+zeros(size(x))'];
    end
    %replace par1 par2 with par(1) par(2) etc.
    f = regexprep(f, 'par([0-9])+', 'par($1)');
    fun = str2func(f);
end

function res = run_derivatives(x, pars, funs)
    res = cell(size(funs));
    for i = 1:numel(funs)
        if isnumeric(funs{i})
            res{i} = zeros(size(x)) + funs{i};
        else
            res{i} = zeros(size(x)) + funs{i}(x, pars);
        end
    end
end


function [cost, LL] = Cost_Parametric(par, J, N, X0, X, dt, EZfun, H, H_coeff, Coeff, mesh_fine, intsigma, mufuns, sigmafuns, sigma)
    %create spline structs + derivatives
    if ~iscell(mufuns)
        mu_vals = mufuns(X0, par);
        sigma_vals = sigmafuns(X0, par);
    else
    mu_vals = run_derivatives(X0, par, mufuns);
    sigma_vals = run_derivatives(X0, par, sigmafuns);
    end
    sigma_meshfine = sigma(mesh_fine, par);
    sigma_X = sigma(X, par);
    %in the EZ functions we now always include dt 
    EZ = EZfun(dt, mu_vals{:}, sigma_vals{:});
    if all(sigma_meshfine > 0)
        VAR = EZ(2, :) - EZ(1, :).^2;
        if min(VAR) > 0
            rho = VAR.^(-1 ./ 2);
            % cpars = num2cell(par);
            Z = intsigma(X0, X, par) ./ sqrt(dt);
            %alternative to the numeric integral of Splines
            %     z = cumtrapz(mesh_fine, 1 ./ sigma_meshfine);
            %     z = @(x)interp1(mesh_fine,z,x);
            %     Z = (z(X) - z(X0)) ./ sqrt(dt); %dt=1 is chosen in the command Z=(z(X)-z(X0))./sqrt(dt);
            Z1 = rho .* (Z - EZ(1, :));
            EZh = zeros(J, N - 1);
            EZh(2, :) = 1;
            for j = 3:J
                EZh(j, :) = rho.^j .* sum([EZ(j: -1:1, :).' ones(N - 1, 1)] .* Coeff{j}(EZ(1, :).'), 2).';
            end
     
            %Hermite expansion coeficients of transition density for arbitrary number of coefficients (J) and arbitrary
            %temporal taylor expansion order (K) for each coefficient
            eta = zeros(J, N - 1);
            for j = 3:2:J
                eta(j, :) = 1 / shared.factorial1(j) .* sum(H_coeff{j} .* EZh(j: -2:1, :));
            end
            for j = 4:2:J
                eta(j, :) = 1 / shared.factorial1(j) .* sum(H_coeff{j} .* [EZh(j: -2:1, :); ones(1, N - 1)]);
            end
            eta = eta(3:end, :); % We do not need the first 2 rows of eta as they are 0
     
            A = rho ./ (sqrt(dt) .* sigma_X) .* normpdf(Z1); %dt added here
            B = [ones(1, N - 1); eta];
            C = [ones(1, N - 1); H(Z1)];
            D = sum(B .* C);
            L = A .* D; % L is the vector of likelihoods
     
            %The following is a 'death penalty' implementation of the constraint of not geting illegitimate likelihoods happening when parameters are rather 'far' from the optimal parameters
            if ~isreal(L)
                cost = realmax;LL = nan;
            elseif min(L) > 0
                cost = -sum(log(L), 'omitnan'); % objective function is the negative of sum of log-likelihoods
                LL = log(L);
            else
                cost = realmax;LL = nan;
            end
        else
            cost = realmax;LL = nan;
        end
    else
        cost = realmax;LL = nan;
    end
end


