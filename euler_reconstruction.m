%euler_reconstruction - reconstruct Langevin equation using the Euler
%method.
%
%usage:
%res=euler_reconstruction(data,dt,'name',value,...)
%data - vector with the equally spaced time series
%dt   - the time interval between two data points
%name-value pairs:
%'lb'  - vector with the lower bounds of all parameters/knots (defaults for
%        spline knots: -10 for mu and 0 for sigma)
%'ub'  - vector with the lower bounds of all parameters/knots (defaults for
%        spline knots: 10)
%'L'   - left boundary of the domain (default min(data))
%'R'   - right boundary for the domain (default max(data))
%'solver' solver for the log-likelihood (default fmincon)
%       fmincon, GWO and PS
%'solveroptions' object or struct with additional options dependent on the 
%      used solver (for instance made with optimoptions)
%'useparallel' use parallel computing (default false)
%'search_agents', number of search agents for GWO
%'maxiter' maximum number of iterations
%
% Spline (cannot be combined with parametric): 
%'nknots' the number of knots for mu and sigma use [n 1] for additive noise (default [8 8])
%'knots' (alternative for nknots) the values of the knots
%'spline' type of splines for mu and sigma (for instance LL for both
%         linear) or 'SCSL' ('SCS' for mu and 'L' for sigma)  default CC
%        'L' = linear interpolation
%        'C' = cubic (order 2) spline
%        'Q' = order 3 spline
%        'P' = pchip 
%        'SCS' = cubic smoothing spline
%        'Approximate' = SCS but uses L for fast fitting
%'transform' transform data before fitting and transform parameters back
%         values: 'none','zscore'
%
%Parametric (cannot be combined with Spline):
%'mu'    parametric function handle for mu
%'sigma' parametric function handle for sigma (empty = additive)
%'npars' (optional) number of parameters default automatic
%'gradient_fun' handle to gradient function (generated with eulergrad)
%
%
function result = euler_reconstruction(data, dt, varargin)
    defaults_param = struct('prev', [], 'prev_range', 0.05, 'lb', [], 'ub', [], 'l', [], 'r', [], 'npars', [], ...
        'mu', [], 'sigma', [], 'gradient_fun', [], 'hessian_fun', [], 'solver', 'fmincon', 'reconst_fraction', [],'solveroptions', [], 'useparallel', false, 'zscore', false, 'search_agents', 1, 'maxiter', realmax);
    defaults_spline = struct('prev', [], 'prev_range', 0.05, 'lb', [], 'ub', [], 'l', [], 'r', [], 'nknots', [8 8],...
        'knots', [], 'spline', 'CC', 'zscore', false, 'gradient_fun', [], 'hessian_fun', [], ...
        'solver', 'fmincon', 'reconst_fraction', [], 'solveroptions', [], 'useparallel', false, 'search_agents', 1, 'maxiter', realmax);

    if ischar(data) && strcmp(data, '-d') %get default options;
        result = struct('parametric', defaults_param, 'spline', defaults_spline);
        return;
    end
    [options, parametric, index_mu, knots] = shared.parse_options_reconst(data, varargin, defaults_param, defaults_spline);

    if options.zscore && parametric
        error('Z-score transformation is not supported for parametric models');
    end
    
    %replace points outside the domain borders (L and R) with NaN
    if ~isempty(options.l) && ~isnan(options.l)
        ndx = data < options.l;
    else
        ndx = false(size(data));
    end
    if ~isempty(options.r) && ~isnan(options.r)
        ndx = ndx | data > options.r;
    end
    data(ndx) = NaN;
 
    if all(isnan(data))
        error('Not enough non-missing data for reconstructions');
    end

    datastats = struct('range', [min(data,[], 'omitnan'), max(data,[], 'omitnan')], ...
        'mean', mean(data, 'omitnan'), 'std', std(data, 'omitnan'));
    if options.zscore
        [data, knots, options] = shared.ztransform(data, datastats, knots, options);
    end

    % Handling big data
    fraction = options.reconst_fraction;
    if isempty(fraction)
        X0 = data(1:end - 1);
        X = data(2:end);
    else
        strata = fraction(1);
        p = fraction(2);
        sample_pairs = sample_stratified(data,strata,p);
        writematrix(sample_pairs,'sample_pairs.csv');
        X0 = sample_pairs(1:2:end);
        X = sample_pairs(2:2:end);
    end

    if size(X0, 1) > 1
        X0 = X0.';
        X = X.';
    end

    hfig = findobj(0, 'tag', 'optimdlg');
    if ~isempty(hfig)
        close(hfig);
    end
    tim = tic;
    %in this loop we try to find initial parameters such that the sigma function 
    %is fully positive. This is needed for fmincon to succeed
    fine_mesh = linspace(min(data, [], 'omitnan'), max(data, [], 'omitnan'), 5000);
    postive_sigma = false;

    if ~parametric
        spline_mu = options.spline{1};
        spline_sigma = options.spline{2};
        iter = 1;
        [cost, mu, sigma] = Create_Euler_Cost_Spline(X0, X, dt, knots, index_mu, spline_mu, spline_sigma);
        while ~postive_sigma && iter < 100
            par0 = options.lb + (options.ub - options.lb) .* rand(1, length(options.lb));
            iter = iter + 1;
            if any(strcmp(options.solver, {'fmincon', 'fminunc'})) && ~parametric && any(strcmp(spline_sigma, {'Q', 'C', 'SCS', 'Approximate'}))
                %these spline types can cause an endless loop in rare cases
                %therefore we limit the initial par0 to values around the
                %optimal fit for the linear interpolation model (LL)
                cost_linear = Create_Euler_Cost_Spline(X0, X, dt, knots, index_mu, 'L', 'L');
                %Create_Euler_Cost_Funcs(X0, X, dt, knots, M, ModelType, 'LL');
                %fast fmincon for initial guess
                par_linear = fmincon(cost_linear, par0, [], [], [], [], options.lb, options.ub, []);
                % we choose a random parameter value between 0.5*parLL and
                % 1.5*parLL
                par0 = par_linear * 0.5 + par_linear .* rand(size(par0));
            end
          
            %Calculating the objective function 
            %cost = @(par)Cost_Spline_Euler(par,X0,X,knots,M,ModelType,SplineType);
            postive_sigma = any(strcmp(spline_sigma, {'L', 'P'})) || all(sigma(fine_mesh, par0) > 0);
        end

    else %parametric
        iter = 1;
        while ~postive_sigma && iter < 100
            par0 = options.lb + (options.ub - options.lb) .* rand(1, length(options.lb));
            iter = iter + 1;
         
            %Calculating the objective function 
            %cost = @(par)Cost_Spline_Euler(par,X0,X,knots,M,ModelType,SplineType);
            postive_sigma = all(options.sigma(fine_mesh, par0) > 0);
        end
        cost = Create_Cost_Funcs_Parameteric(X0, X, dt, options.mu, options.sigma, options.gradient_fun);
    end

    dim = length(par0);
    switch lower(options.solver)
        case 'fmincon'
            fprintf('\n Solving maximum likelihood problem (Euler)\n');
            %Optimization using fmincon
            if parametric && isempty(options.gradient_fun)
                opts = optimoptions('fmincon', 'FiniteDifferenceStepSize', 10^(-6), 'UseParallel', false, 'MaxFunEvals', realmax, 'MaxIter', options.maxiter, ...
                    'OutputFcn', @shared.optimdlg);
            end
            if parametric && ~isempty(options.gradient_fun)
                opts = optimoptions('fmincon', 'FiniteDifferenceStepSize', 10^(-6), 'UseParallel', options.useparallel, 'MaxFunEvals', realmax, 'MaxIter', options.maxiter, ...
                    'SpecifyObjectiveGradient', parametric && ~isempty(options.gradient_fun),'OutputFcn', @shared.optimdlg);
            end

            if ~parametric && isempty(options.gradient_fun)
                opts = optimoptions('fmincon', 'FiniteDifferenceStepSize', 10^(-6), 'UseParallel', false, 'MaxFunEvals', realmax, 'MaxIter', options.maxiter, ...
                    'OutputFcn', @shared.optimdlg);
            end

            lb=options.lb;ub=options.ub;
            opts = optimoptions(opts,'display','none');

            if options.search_agents == 1
                [estimated_par, f_best] = fmincon(cost, par0, [], [], [], [], options.lb, options.ub, [], opts);
            else
                opts=optimoptions(opts,'display','none');
                p=sobolset(dim,'Skip',1e3,'Leap',1e2);p=scramble(p,'MatousekAffineOwen');  % sobol random sampling (is the best especially at high dimensions)
                RS=net(p,options.search_agents);   % This is a random sample whose entities are in (0,1)

                Estimated_par=zeros(options.search_agents,dim);F_best=zeros(options.search_agents,1);
                A=[];
                while size(A,1)<options.search_agents   % this is to make sure we meet the requirements of user to get >= search_agents legitimate solutions
                    ptmatrix=lb+(ub-lb).*RS;
                    parfor n=1:options.search_agents
                        [Estimated_par(n,:), F_best(n)] = fmincon(cost, ptmatrix(n,:), [], [], [], [], lb, ub, [], opts);
                    end
                    F_best=F_best(:);
                    A=[A;Estimated_par F_best];
                    idx=A(:,end)>1.7900e+308;
                    A=A(~idx,:);
                end
                A=sortrows(A,dim+1);
                estimated_par=A(1,1:dim);f_best=A(1,end);
            end
        case 'gwo'
            fprintf('\n Solving maximum likelihood problem (Euler)\n');
            %Optimization using Improved Grey-Wolf optimizer
            if parametric
                [~, mundx] = shared.getfundimension(options.mu, []);
                M = sum(mundx);
                %find(arrayfun(@(x)has(sigma,par(x)), 1:length(par)) == 1, 1);
                [estimated_par, f_best] = shared.IGWO_Parametric(dim, options.search_agents, options.maxiter, options.lb, options.ub, cost, M, options.useparallel, @shared.optimdlg);
            else
                M = options.nknots(1);
                switch options.maxiter
                    case realmax
                        % if strcmpi(graph, 'Yes')
                        %     [estimated_par, f_best] = IGWO_Spline(dim, M, SearchAgents_no, options.maxiter, options.lb, options.ub, cost, dt, options.useparallel, graph, knots_mu, mesh_fine, ModelType, SplineType); %Results (of each iteration) are displayed and ploted
                        % else
                        [estimated_par, f_best] = shared.IGWO_Spline_Euler(dim, M, options.search_agents, options.maxiter, options.lb, options.ub, cost, options.useparallel, @shared.optimdlg); %Results are displayed only
                        % end
                    otherwise %Results are neighter displayed nor ploted
                        [estimated_par, f_best] = shared.IGWO_Spline_Euler(dim, M, options.search_agents, options.maxiter, options.lb, options.ub, cost, options.useparallel, @shared.optimdlg);
                        if isequal(f_best, realmax)
                            disp('The algorithm did not find a legitimate solution. Increase maximum iteration number ''maxiter'' and if this does not work consider simpler models with less parameters');
                        end
                end
            end
        otherwise
            error('Solver ''%s'' is not supported, valid options are ''fmincon'',''GWO'',''PS''', options.solver);
    end
    if isequal(f_best, realmax)
        disp('The algorithm did not find a legitimate solution. Increase maximum iteration number options.maxiter and if this does not work consider simpler models with less parameters');
    else
        disp('Estimated parameters : ');
        disp(num2str(estimated_par));
        disp(['- sum of log-likelihoods) : ' num2str(f_best)]);
    end
    hfig = findobj(0, 'tag', 'optimdlg');
    if ~isempty(hfig)
        close(hfig);
    end
    %   fine_mesh = linspace(min(data), max(data), 5000);
    tim2 = toc(tim);
    if ~parametric
        if any(sigma(fine_mesh, estimated_par) < 0)
            error('The fitted sigma function has negative points, please refit or use higher lower boundaries');
        end

        if strcmp(spline_sigma, 'Approximate')
            costSCS = Create_Euler_Cost_Spline(X0, X, dt, knots, index_mu, spline_mu, 'SCS');
            f_best = costSCS(estimated_par);
            fprintf('f_best (real)= %g\n', f_best);
        end

        result = struct('estimated_par', estimated_par, 'negloglik', f_best, 'datarange', [min(data), max(data)], 'dt', dt, ...
            'mufun', mu, 'sigmafun', sigma, 'costfun', cost, 'computation_time', tim2, 'knots', {{knots(index_mu), knots(~index_mu)}}, ...
            'options', options);
        if options.zscore
            result = shared.back_ztransform(result);
        end
    else
        %paramtric has empty knots
        if ~strcmp(options.solver, 'fitnlm')
            STD1 = options.sigma(fine_mesh, estimated_par);
            if any(STD1 < 0)
                error('The fitted sigma function has negative points, please refit or use higher lower boundaries');
            end
            if numel(STD1) == 1
                %if the sigma function returns a single value we adapt it to
                %return a vector of the same size as x
                options.sigma = @(x,par)zeros(size(x))+options.sigma(x,par);
            end
        end
        result = struct('estimated_par', estimated_par, 'negloglik', f_best, 'datarange', [min(data), max(data)], 'dt', dt, ...
            'mufun', options.mu, 'sigmafun', options.sigma, 'computation_time', tim2, 'knots', [], ...
            'options', options);
    end

end

function costfun = Create_Cost_Funcs_Parameteric(X0, X, dt, mu, sigma, gradient_fun)
    function [cost, g_cost] = cost_with_grad(par, X0, X, dt, sqrtdt, mu, sigma, gradient_fun)
        STD1 = sigma(X0, par);
        cost = -sum(-0.5 .* (log(2 .* pi .* dt .* STD1.^2) + ((X - X0 - mu(X0, par) .* dt) ./ (STD1 .* sqrtdt)).^2), 'omitnan');
        if nargout>1
            cpar = num2cell(par);
            g_cost = -sum(gradient_fun(X0, X, dt, cpar{:}),2,'omitnan');
        end
    end
    function cost = cost_wo_grad(par, X0, X, dt, sqrtdt, mu, sigma)
        STD1 = sigma(X0, par); %done once for efficiency
        cost = -sum(-0.5 * (log(2 .* pi .* dt .* STD1.^2) + ((X - X0 - mu(X0, par) .* dt) ./ (STD1 .* sqrtdt)).^2), 'omitnan');
    end
    sqrtdt = sqrt(dt);
    if isempty(gradient_fun)
        costfun = @(par)cost_wo_grad(par,X0,X,dt,sqrtdt,mu,sigma);
    else
        costfun = @(par)cost_with_grad(par,X0,X,dt,sqrtdt,mu,sigma,gradient_fun);
    end
end

function cost = euler_likelihood(X, X0, mu, sigma, dt, sqrtdt)
    %this is the base formula to calculate the -log(likelihood) using the
    %Euler method
    %normal distibution with mean X0-mu*dt and standard deviation
    %sqrt(dt)*sigma
    % 
        mu1=X0+mu.*dt;
        sigma1=sigma.*sqrtdt;
        z = (X - mu1) ./ sigma1;
        L = -.5 .*z.*z - log(sqrt(2 .*pi).*sigma1);
        cost=-sum(L,'omitnan');

    %slightly less operations than
    % cost = -sum(-0.5 .* (log(2 .* pi .* dt .* sigma.^2) + ((X - X0 - mu .* dt) ./ (sigma .* sqrtdt)).^2), 'omitnan');
end

function [costfun, mu, sigma] = Create_Euler_Cost_Spline(X0, X, dt, knots, index_mu, spline_mu, spline_sigma)
    % we create an efficient costfunction and the used mu and sd functions
    function cost = ffun(par, X0, X, dt, sqrtdt, sigma, mu)
        %much faster than the functions with fnmin, I assume that if all sigma's are positive, fnmin is
        %also positive
        STD1 = sigma(X0, par);
        if all(isnan(STD1) | STD1 > 0)
            cost = euler_likelihood(X, X0, mu(X0, par), STD1, dt, sqrtdt);
        else
            cost = realmax;
        end
    end
    function cost = ffun_approx(par, X0, X, dt, sqrtdt, sigma, mu, testsigma)
        %Version for the approximate type, it uses another function for
        %testing positivity
        if all(isnan(X0) | testsigma(X0, par) > 0)
            cost = euler_likelihood(X, X0, mu(X0, par), sigma(X0, par), dt, sqrtdt);
        else
            cost = realmax;
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
    knots_mu = knots(index_mu);
    knots_sigma = knots(index_sigma);
    sqrtdt = sqrt(dt);
    %the indices to mu and sigma parameters   
 
    %the parameters for mu
    switch spline_mu
        case 'L'
            mu = @(x,par)interp1(knots_mu,par(index_mu),x);
        case 'Q'
            mu = @(x,par)fnval(spapi(3, knots_mu, par(index_mu)),x);
        case 'C'
            mu = @(x,par)ppval(spline(knots_mu, par(index_mu)),x);
            %  mu = @(x,par)interp1(knots_mu,par(ndx),x,'cubic');
        case 'P'
            %mu = @(x,par)interp1(knots_mu,par(ndx),x,'pchip');
            mu = @(x,par)ppval(pchip(knots_mu, par(index_mu)),x);
        case {'SCS', 'Approximate'}
            mu = @(x,par)ppval(csaps(knots_mu, par(index_mu)),x);
        case 'constant'
            mu = @(x,par)par(index_mu)+zeros(size(x));
        otherwise
            error('Spline type of mu ''%s'' is not supported, valid options are ''Approximate'',''L'',''Q'',''C'',''P'',''SCS''', spline_mu)
    end

    %the parameters for sigma
    switch spline_sigma
        case 'constant'
            %additive noise
            sigma = @(x,par)par(index_sigma)+zeros(size(x));
            costfun = @(par)euler_likelihood(X, X0, mu(X0, par), par(index_sigma), dt, sqrtdt);
        case 'L'
            sigma = @(x,par)interp1(knots_sigma,par(index_sigma),x);
            costfun = @(par)euler_likelihood(X, X0, mu(X0, par), sigma(X0, par), dt, sqrtdt);
        case 'Q'
            sigma = @(x,par)fnval(spapi(3, knots_sigma, par(index_sigma)),x);
            costfun = @(par)ffun(par, X0, X,dt,sqrtdt, sigma, mu);
        case 'C'
            sigma = @(x,par)ppval(spline(knots_sigma, par(index_sigma)),x);
            costfun = @(par)ffun(par, X0, X,dt,sqrtdt, sigma, mu);
            %sigma = @(x,par)interp1(knots_sigma,par(1:M),x,'cubic');
        case 'P'
            sigma = @(x,par)ppval(pchip(knots_sigma, par(index_sigma)),x);
            costfun = @(par)euler_likelihood(X, X0, mu(X0, par), sigma(X0, par), dt, sqrtdt);
        case 'Approximate'
            sigma = @(x,par)ppval(csaps(knots_sigma, par(index_sigma)),x);
            sigma_linear = @(x,par)interp1(knots_sigma,par(index_sigma),x);
            mu_linear = @(x,par)interp1(knots_mu,par(index_mu),x);
            costfun = @(par)ffun_approx(par, X0, X,dt,sqrtdt, sigma_linear, mu_linear,sigma);
        case 'SCS'
            sigma = @(x,par)ppval(csaps(knots_sigma, par(index_sigma)),x);
            costfun = @(par)ffun(par, X0, X, dt,sqrtdt,sigma, mu);
        otherwise
            error('Spline type of sigma ''%s'' is not supported, valid options are ''Approximate'',''L'',''Q'',''C'',''P'',''SCS''', spline_sigma)
    end

end



