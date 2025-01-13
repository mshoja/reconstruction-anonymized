function [options, parametric, index_mu, knots] = parse_options_reconst(data, args, defaults_param, defaults_spline)
    %this function parses the options for both euler_reconprestruction and
    %hermite_reconstruction (as this largely overlaps)
    opts = reshape(args, 2, length(args) / 2);
    opts(1, :) = lower(opts(1, :)); %case insensitive
    ndx = strcmp('prev', opts(1, :));
    if any(ndx)
        prev_results = opts{2, ndx};

    else
        prev_results = [];
    end
    knots=[];
    ndx = strcmp('mu', opts(1, :));
    parametric = (any(ndx) && ~isempty(opts{2, ndx})) || (~isempty(prev_results) && isfield(prev_results.options, 'mu'));
    %default options
    if parametric
        options = defaults_param;
        if any(ismember(opts(1, :), {'spline', 'knots', 'nknots'}))
            error('You cannot combine spline parameters (spline,knots,nknots) with paramatric (mu,sigma,npars)');
        end
    else
        options = defaults_spline;
        if any(ismember(opts(1, :), {'mu', 'sigma', 'npars'}))
            error('You cannot combine spline parameters (spline,knots,nknots) with paramatric (mu,sigma,npars)');
        end
    end
    validoptions = fieldnames(options);
    %give error if there are unknown options
    if any(~ismember(opts(1, :), validoptions))
        s1 = sprintf('"%s",', opts{1, ~ismember(opts(1, :), validoptions)});
        s2 = sprintf('%s, ', validoptions{:});
        error('fitlangevin:options', 'Unknown options in reconstruction (%s)\nValid options are: %s\n', s1(1:end - 1), s2);
    end
    %change the default settings for the number of "search agents"
    ndx = strcmp(opts(1, :), 'solver');
    if ~any(ismember(opts(2, ndx), {'fmincon', 'fminunc'}))
        options.search_agents = 5;
    end
    %update the options struct
    if ~isempty(prev_results)
        %replace the default settings with the options used for the
        %euler_reconstruction
        opt_prev = prev_results.options;
        f = fieldnames(opt_prev);
        for i = 1:numel(f)
            %only if the field name is valid we set the option
            ndx = strcmp(f{i}, validoptions);
            if any(ndx)
                options.(f{i}) = opt_prev.(f{i});
            end
        end
        %do not use the previous ub and lb
        if ~isempty(prev_results.estimated_par)
           options.ub = [];
           options.lb = [];
        end
    end
    %override the previous results with the options supplied by the user
    for i = 1:numel(validoptions)
        opt1 = validoptions{i};
        ndx = strcmp(opts(1, :), opt1);
        if any(ndx)
            options.(opt1) = opts{2, ndx};
        end
    end
    if ~parametric
        if ~isempty(options.nknots)
            if numel(options.nknots) == 1
                options.nknots = [options.nknots options.nknots];
            end
            mmin = min(data,[],'omitnan');
            if isfield(options,'l')&&~isempty(options.l)&&~isnan(options.l)
                mmin=options.l;
            end
            mmax = max(data,[],'omitnan');
            if isfield(options,'r')&&~isempty(options.r)&&~isnan(options.r)
                mmax=options.r;
            end
            knots_mu = linspace(mmin, mmax, options.nknots(1));
            knots_sigma = linspace(mmin, mmax, options.nknots(2));
        end
        if iscell(options.knots)
            knots_mu = options.knots{1};
            if numel(options.knots) > 1
                knots_sigma = unique(options.knots{2});
            else
                knots_sigma = unique(options.knots{1});
            end
            options.knots = {knots_mu knots_sigma};
        elseif ~isempty(options.knots)
            knots_mu = unique(options.knots);
            knots_sigma = unique(options.knots);
            options.knots = {knots_mu knots_sigma};
        end
        if ischar(options.spline)
            options.spline = regexp(options.spline, '(Approximate|SCS|C|L|P|Q||[a-zA-Z]*)', 'match');
            if length(options.spline) > 2
                error('cannot read the spline types "%s"', sprintf('%s', options.spline{:}))
            end
        end
        if length(options.spline) == 1
            options.spline = repmat(options.spline, 1, 2);
        end
        if isempty(options.spline)
            options.spline = {[], []};
        end
        M = length(knots_mu);
        A = options.lb(M + 1:end);
        A(A <= 1E-10) = 1E-10;
        options.lb(M + 1:end) = A; %In some rare cases 'eps' does not work, so we chose 10^(-10)
        npars = numel(knots_mu) + numel(knots_sigma);
        index_mu=[true(1,numel(knots_mu)) false(1,numel(knots_sigma))];
        knots=[knots_mu knots_sigma];
    else
        [npars,index_mu] = shared.getfundimension(options.mu, options.sigma);
        if ~isempty(options.npars)
            npars = options.npars;            
        end
        if isempty(options.sigma)||ischar(options.sigma)&&strcmpi(options.sigma,'additive')
            %additive noise by default
            options.sigma = @(x,par)par(npars)+zeros(size(x));
        end
    end
%  
%     if isempty(prev_results)
%         options.prev_range = [];
%     end
    if ~isempty(prev_results) && ~isempty(prev_results.estimated_par) && numel(prev_results.estimated_par) == npars

        if parametric
            % prev_lb = prev_results.estimated_par - prev_results.estimated_par * options.prev_range;
            % prev_ub = prev_results.estimated_par + prev_results.estimated_par * options.prev_range;

            ran = (max(prev_results.estimated_par) - min(prev_results.estimated_par)) / 2;
            prev_lb = prev_results.estimated_par - ran;
            prev_ub = prev_results.estimated_par + ran;
           
        else
            %if bounds are not set by the user, use range close to the prev results
            %spline version
            % if isempty(options.prev_range)
            %     options.prev_range=0.05;
            % end
            % prev_lb = prev_results.estimated_par;
            % prev_ub = prev_lb;
            % 
            % ind = 1:npars;
            % ndx_sigma = ind > numel(knots_mu);
            % ran = (max(prev_results.estimated_par(~ndx_sigma)) - min(prev_results.estimated_par(~ndx_sigma))) / 2;
            % 
            % if ran == 0
            %     ran = max(prev_results.estimated_par(~ndx_sigma)) / 2;
            % end
            % prev_lb(~ndx_sigma) = prev_results.estimated_par(~ndx_sigma) - ran * options.prev_range;
            % prev_ub(~ndx_sigma) = prev_results.estimated_par(~ndx_sigma) + ran * options.prev_range;
            % ran = (max(prev_results.estimated_par(ndx_sigma)) - min(prev_results.estimated_par(ndx_sigma))) / 2;
            % 
            % if ran == 0
            %     ran = max(prev_results.estimated_par(ndx_sigma)) / 2;
            % end
            % prev_lb(ndx_sigma) = prev_results.estimated_par(ndx_sigma) - ran * options.prev_range;
            % prev_ub(ndx_sigma) = prev_results.estimated_par(ndx_sigma) + ran * options.prev_range;
            % 
            % ndx_sigma1 = ndx_sigma & prev_lb < eps;
            % prev_lb(ndx_sigma1) = eps;

            % prev_lb = prev_results.estimated_par;
            % prev_ub = prev_lb;

            ran = (max(prev_results.estimated_par) - min(prev_results.estimated_par)) / 2;
            prev_lb = prev_results.estimated_par - ran;
            prev_ub = prev_results.estimated_par + ran;

            ind = 1:npars;
            ndx_sigma = ind > numel(knots_mu);
            ndx_sigma1 = ndx_sigma & prev_lb < eps;
            prev_lb(ndx_sigma1) = eps;
        end
    end
    if isempty(options.lb)
        if ~parametric
            if ~isempty(prev_results) && numel(prev_results.estimated_par) == npars
                %if bounds are not set by the user, use range close to the prev results
                %spline version 
                options.lb = prev_lb;
            else
                options.lb = -10 + zeros(1, npars);
                options.lb(numel(knots_mu) + 1:end) = eps;
            end
        else
            if ~isempty(prev_results)
                options.lb = prev_lb;
            else
                options.lb = -100 + zeros(1, npars);
            end
        end
    end
    if isempty(options.ub)
        if ~parametric
            if ~isempty(prev_results)
                options.ub = prev_ub;
            else
                options.ub = 10 + zeros(1, npars);
            end
        else
            if ~isempty(prev_results)
                options.ub = prev_ub;
            else
                options.ub = 100 + zeros(1, npars);
            end
        end
    end
    
    if numel(options.lb) ~= npars || numel(options.ub) ~= npars
        error('The size of the lower bound and upper bound should be the same as the number of parameters %d', npars);
    end
end
