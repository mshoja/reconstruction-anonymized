function [costfun, mu1, sigma1, EZ] = Create_Hermite_Cost_Spline(data, dt, J, K, knots, index_mu, spline_mu, spline_sigma, refined, fraction)
    % we create an efficient costfunction and the used mu and sd functions
    function res = normal_fun(aspline_derivs, x)
        s = aspline_derivs;
        res = spline_val(s{1}, x);
    end

    function nder = num_derivatives(aspline)
        switch aspline
            case {'L', 'Approximate'}
                nder = 1;
            case 'Q'
                nder = 2;
            case {'C', 'P', 'SCS'}
                nder = 3;
            case 'constant'
                nder = 0;
            otherwise
                nder = 3;
        end
    end

    if sum(index_mu) == 1
        spline_mu = 'constant';
        %same mu everywhere 
    end
    index_sigma = ~index_mu;
    if sum(index_sigma) == 1
        spline_sigma = 'constant';
        %additive noise 
    end

    knots_mu=knots(index_mu);
    knots_sigma=knots(index_sigma);
    
    nder_mu = num_derivatives(spline_mu);
    nder_sigma = num_derivatives(spline_sigma);
    %if we have a different number of derivatives for mu and sigma and we
    %have not additive noise (=nder_sigma>0) then we use the maximum of
    %both. This saves creating a lot of different EZ functions if we
    %combine different kinds of splines for mu and sigma. 
    if nder_mu ~= nder_sigma && nder_sigma > 0
        nder_mu = max(nder_mu, nder_sigma);
        nder_sigma = nder_mu;
    end
    %the parameters for mu
    switch spline_mu
        case {'L', 'Approximate'}
            mu_spline = @(par)spline_derivs(interp1(knots_mu, par(index_mu), 'linear', 'pp'),nder_mu); %#ok<INTRPP>
        case 'Q'
            mu_spline = @(par)spline_derivs(spapi(3, knots_mu, par(index_mu)),nder_mu);
        case 'C'
            mu_spline = @(par)spline_derivs(spline(knots_mu, par(index_mu)),nder_mu);
        case 'P'
            mu_spline = @(par)spline_derivs(pchip(knots_mu, par(index_mu)),nder_mu);
        case 'SCS'
            mu_spline = @(par)spline_derivs(csaps(knots_mu, par(index_mu)),nder_mu);
        case 'constant'
            mu_spline = @(par)spline_derivs(par(index_mu),nder_mu);
        otherwise
            error('Spline type of mu ''%s'' is not supported, valid options are ''Approximate'',''L'',''Q'',''C'',''P'',''SCS''', spline_mu)
    end
    
    %the parameters for sigma
    switch spline_sigma
        case 'constant'
            %additive noise
            sigma_spline = @(par)spline_derivs(par(index_sigma),nder_sigma);
        case 'L'
            sigma_spline = @(par)spline_derivs(interp1(knots_sigma, par(index_sigma), 'linear', 'pp'),nder_sigma); %#ok<INTRPP>
        case 'Q'
            sigma_spline = @(par)spline_derivs(spapi(3, knots_sigma, par(index_sigma)),nder_sigma);
        case 'C'
            sigma_spline = @(par)spline_derivs(spline(knots_sigma, par(index_sigma)),nder_sigma);
        case 'P'
            sigma_spline = @(par)spline_derivs(pchip(knots_sigma, par(index_sigma)),nder_sigma);
        case 'Approximate'
            sigma_spline = @(par)spline_derivs(interp1(knots_sigma, par(index_sigma), 'linear', 'pp'),nder_sigma);
        case 'SCS'
            sigma_spline = @(par)spline_derivs(csaps(knots_sigma, par(index_sigma)),nder_sigma);
        otherwise
            error('Spline type of sigma ''%s'' is not supported, valid options are ''Approximate'',''L'',''Q'',''C'',''P'',''SCS''', spline_sigma)
    end
    mu1 = @(x,par)normal_fun(mu_spline(par),x);
    sigma1 = @(x,par)normal_fun(sigma_spline(par),x);
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

    if refined %refined Aït-Sahalia
        %mesh_fine is used to approximate the integral of 1/sigma
        mesh_fine = linspace(min(data), max(data), 2000);
        %we do not need the first two elements
        [H, H_coeff] = shared.Hfun(J, true);
        costfun = @(par)Cost_Spline_refined(par, J, N, X0, X, dt, EZ, H, H_coeff, shared.Coefffun(J),mesh_fine ,mu_spline,sigma_spline);
    else %original Aït-Sahalia
        H = shared.Hfun(J, false);
        costfun = @(par)Cost_Spline_AitSahalia(par, N, X0, X, dt, EZ, H, mu_spline, sigma_spline);
    end
end

function cost = Cost_Spline_AitSahalia(par, N, X0, X, dt, ETAfun, H, mu_spline, sigma_spline)
    %original Aït-Sahalia method for splines, note that we can only use additive noise
    %in this case
    mu_deriv = mu_spline(par);
    mu_vals = cell(size(mu_deriv));
    for i = 1:length(mu_deriv)
        %run the splines for all derivatives
        mu_vals{i} = spline_val(mu_deriv{i}, X0);
    end
    sigma_deriv = sigma_spline(par);
    sigma_vals = cell(size(sigma_deriv));
    for i = 1:length(sigma_deriv)
        sigma_vals{i} = spline_val(sigma_deriv{i}, X0);
    end
    if length(sigma_deriv) > 1
        error('For splines, the Aït-Sahalia method only can deal with additive noise, set ''refined'' true to fit multiplicative noise')
    end
    sigma = mean(sigma_vals{1});
    s = 1 / sigma;
    for i = 1:length(mu_deriv)
        mu_vals{i} = mu_vals{i} .* s;
        s = s * sigma;
    end
    
    %in the ETA functions we now always include dt
    ETA = ETAfun(dt, mu_vals{:});
    sqrt_dt_sigma=sigma.*sqrt(dt);
    Z = (X - X0) ./ sqrt_dt_sigma; %dt=1 is chosen in the command diff_data./(p_end*sqrt(dt));
    A = 1 ./ sqrt_dt_sigma.* normpdf(Z); %dt=1 is chosen in the command A=1 ./(sqrt(dt).*p_end).*normpdf(Z);
    B = [ones(1, N - 1); ETA];
    C = [ones(1, N - 1); H(Z)];
    D = sum(B .* C);
    L = A .* D; % L is the vector of likelihoods
    L = L(~isnan(L)); %this is to account for replicate and missing data

    %The following is a 'death penalty' implementation of the constraint of not geting illegitimate likelihoods happening when parameters are rather 'far' from the optimal parameters
    if ~isreal(L)
        cost = realmax;
    elseif min(L) > 0
        cost = -sum(log(L), 'omitnan');   % objective function is the negative of sum of log-likelihoods
    else
        cost = realmax;
    end
end


function cost = Cost_Spline_refined(par, J, N, X0, X, dt, EZfun, H, H_coeff, Coeff, mesh_fine, mu_spline, sigma_spline)
    %refined Aït-Sahalia method. Can also be used for multiplicative noise
    %create spline structs + derivatives
    mu_deriv = mu_spline(par);
    mu_vals = cell(size(mu_deriv));
    for i = 1:length(mu_deriv)
        %run the splines for all derivatives
        mu_vals{i} = spline_val(mu_deriv{i}, X0);
    end
    sigma_deriv = sigma_spline(par);
    sigma_vals = cell(size(sigma_deriv));
    for i = 1:length(sigma_deriv)
        sigma_vals{i} = spline_val(sigma_deriv{i}, X0);
    end
    sigma_meshfine = spline_val(sigma_deriv{1}, mesh_fine);
    sigma_X = spline_val(sigma_deriv{1}, X);
    %in the EZ functions we now always include dt 
    EZ = EZfun(dt, mu_vals{:}, sigma_vals{:});
    if all(sigma_meshfine > 0)
        VAR = EZ(2, :) - EZ(1, :).^2;
        if min(VAR) > 0
            rho = VAR.^(-1 ./ 2);
            %alternative to the analytic integral of Parametric ("complex
            %noise")
            %Z = int(1 / sigma, x);
            %Z = (Z - subs(Z, x, x0)) / sqrt(dt);
            %Z = simplify(Z);
            z = cumtrapz(mesh_fine, 1 ./ sigma_meshfine);
            z = @(x)interp1(mesh_fine,z,x);
            Z = (z(X) - z(X0)) ./ sqrt(dt); %dt=1 is chosen in the command Z=(z(X)-z(X0))./sqrt(dt);
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
                cost = realmax;
            elseif min(L) > 0
                cost = -sum(log(L), 'omitnan'); % objective function is the negative of sum of log-likelihoods
            else
                cost = realmax;
            end
        else
            cost = realmax;
        end
    else
        cost = realmax;
    end
end

function res = spline_val(aspline, x)
    if ~isstruct(aspline)
        res = zeros(size(x)) + aspline;
    elseif strcmp(aspline.form, 'pp')
        res = ppval(aspline, x);
    else
        res = fnval(aspline, x);
    end
end

function res = spline_derivs(aspline, nder)
    res = cell(nder + 1, 1);
    res{1} = aspline;
    if isstruct(aspline)
        for i = 1:nder
            aspline = fnder(aspline);
            if all(aspline.coefs == 0)
                res{i + 1} = 0;
            else
                res{i + 1} = aspline;
            end
        end
    else
        %constant
        for i = 1:nder
            res{i + 1} = 0;
        end
    end
end
