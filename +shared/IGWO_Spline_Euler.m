%___________________________________________________________________%
%  An Improved Grey Wolf Optimizer for Solving Engineering          %
%  Problems (I-GWO) source codes version 1.0                        %
%                                                                   %
%  Developed in MATLAB R2018a                                       %
%                                                                   %
%  Author and programmer: M. H. Nadimi-Shahraki, S. Taghian, S. Mirjalili                                           %
%                                                                   %
% e-Mail: nadimi@ieee.org, shokooh.taghian94@gmail.com,   ali.mirjalili@gmail.com                                   %
%                                                                   %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                                                  %
%                                                                   %
%   Main paper: M. H. Nadimi-Shahraki, S. Taghian, S. Mirjalili     %
%               An Improved Grey Wolf Optimizer for Solving         %
%               Engineering Problems , Expert Systems with          %
%               Applicationsins, in press,                          %
%               DOI: 10.1016/j.eswa.2020.113917                     %
%___________________________________________________________________%
%___________________________________________________________________%
%  Grey Wold Optimizer (GWO) source codes version 1.0               %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, A. Lewis             %
%               Grey Wolf Optimizer, Advances in Engineering        %
%               Software , in press,                                %
%               DOI: 10.1016/j.advengsoft.2013.12.007               %
%                                                                   %
%___________________________________________________________________%

% Improved Grey Wolf Optimizer (I-GWO)
function [par, Alpha_score] = IGWO_Spline_Euler(varargin)

    %*****************************************
    %My additions
    dim = varargin{1};
    M = varargin{2};
    N = varargin{3};
    Max_iter = varargin{4};
    lb = varargin{5};
    ub = varargin{6};
    fobj = varargin{7};
    %dt = varargin{8};
    UseParallel = varargin{8};
    if nargin == 9
        outputFcn = varargin{9};
    else
        outputFcn = [];
    end
    %if nargin > 9
    %    knots = varargin{10};
    %    mesh_fine = varargin{11};
    %    ModelType = varargin{12};
    %    SplineType = varargin{13};
    %end

    if UseParallel
        numWorkers = Inf;
    else
        numWorkers = 1;
    end
    %****************************************

    % lu = [lb .* ones(1, dim); ub .* ones(1, dim)];

    % Initialize alpha, beta, and delta positions
    Alpha_pos = zeros(1, dim);
    Alpha_score = inf; %change this to -inf for maximization problems

    Beta_pos = zeros(1, dim);
    Beta_score = inf; %change this to -inf for maximization problems

    Delta_pos = zeros(1, dim);
    Delta_score = inf; %change this to -inf for maximization problems

    % Initialize the positions of wolves
    % Positions = boundConstraint (Positions, Positions, lu);
    % Positions=initialization(N,dim,ub,lb);

    %***********************************************************
    %My additions
    optimValues = struct('bestfval', NaN, 'bestx', zeros(size(lb)), 'iteration', 0);
    if ~isempty(outputFcn)
        outputFcn(optimValues, 'init');
    end
    % Initialization based on various sampling schemes (this rarely gives legitimate initializations where legitimate means having 'defined objective value')
    % It seems that Sobol sampling is the best especially in high dimensions. Here, dimension is the number of parameters to estimate

    % Sobol sampling (default)
    p=sobolset(dim,'Skip',1e3,'Leap',1e2);p=scramble(p,'MatousekAffineOwen');
    RS=net(p,N);   % This is a random sample whose entities are in (0,1)
    LB=repmat(lb,[N 1]);
    UB=repmat(ub,[N 1]);
    Positions=LB+(UB-LB).*RS;

    % Calculate objective function for each wolf
    Fit = zeros(1, N);
    for i = 1:N
        Fit(i) = fobj(Positions(i, :));
    end

    % Personal best fitness and position obtained by each wolf
    pBestScore = Fit;
    pBest = Positions;

    neighbor = zeros(N, N);
    % Convergence_curve=zeros(1,Max_iter);
    iter = 0; % Loop counter
    %% Main loop
    Fit_GWO = zeros(1, N);
    Fit_DLH = zeros(1, N);
    par = zeros(1, dim); %My additions
    stop = false;
    while iter < Max_iter && ~stop
        for i = 1:N
            fitness = Fit(i);
            % Update Alpha, Beta, and Delta
            if fitness < Alpha_score
                Alpha_score = fitness; % Update alpha
                Alpha_pos = Positions(i, :);

                %***********************************
                %My additions
                par(1:M) = Alpha_pos(1:M);
                par(M + 1:end) = Alpha_pos(M + 1:end); %Converting the estimated parameters to the actual unit of data
                par(M + 1:end) = abs(par(M + 1:end)); %In Euler reconstruction the objective function is invariant under negation of diffusion vector: any legitimate (i.e., positive) vector of diffusion function corresponds with an illegitimate (i.e., negative) solution whose absolute vale is the same
                if ~isempty(outputFcn)
                    optimValues = struct('bestfval', Alpha_score, 'bestx', par, 'iteration', iter);
                    if ~isempty(outputFcn)
                        stop = outputFcn(optimValues, 'iter');
                    else
                        if isequal(Max_iter, realmax) %Results (of each iteration) are displayed only
                            disp(['Iteration = ', num2str( iter)]);
                            disp('Estimated parameters : ');
                            disp(num2str(par));
                            disp(['Approximate value of objective function (negative sum of log-likelihoods): ' num2str(Alpha_score)]);
                            writematrix([par Alpha_score], 'Results.csv', 'WriteMode', 'append');
                        end

                        if nargin > 10 && isequal(Max_iter, realmax) %Results (of each iteration) are displayed and ploted
                            disp(['Iteration = ', num2str( iter)]);
                            disp('Estimated parameters : ');
                            disp(num2str(par));
                            disp(['Approximate value of objective function (negative sum of log-likelihoods): ' num2str(Alpha_score)]);

                            % Grephical illustration for each successful iteration
                            close all;
                            switch ModelType
                                case {'Additive noise1', 'Additive noise2', 'Additive noise'}
                                    switch SplineType
                                        case 'L'
                                            mu = @(x)interp1(knots,par(1:M),x);
                                            plot(mesh_fine, mu(mesh_fine), '-k');
                                            hold on;
                                            plot(mesh_fine, 0 .* mesh_fine, '-k');
                                            hold on;
                                            plot(knots, par(1:M), '*r');
                                        case 'Q'
                                            sp_mu = spapi(3, knots, par(1:M));
                                            plot(mesh_fine, fnval(sp_mu, mesh_fine), '-k');
                                            hold on;
                                            plot(mesh_fine, 0 .* mesh_fine, '-k');
                                            hold on;
                                            plot(knots, par(1:M), '*r');
                                        case {'C', 'P', 'SCS', 'Approximate'}
                                            switch SplineType
                                                case 'C'
                                                    pp_mu = spline(knots, par(1:M));
                                                case 'P'
                                                    pp_mu = pchip(knots, par(1:M));
                                                case {'SCS', 'Approximate'}
                                                    pp_mu = csaps(knots, par(1:M));
                                            end
                                            plot(mesh_fine, ppval(pp_mu, mesh_fine), '-k');
                                            hold on;
                                            plot(mesh_fine, 0 .* mesh_fine, '-k');
                                            hold on;
                                            plot(knots, par(1:M), '*r');
                                    end
                                    annotation('textbox', [0.55, 0.8, 0.1, 0.1], 'String', strcat(['Noise intensity is : ' num2str(par(end))]));
                                    xlabel('State, x');
                                    ylabel('Drift function, \mu (x)');
                                    xlim([-inf inf]);
                                    ylim([-inf inf]);
                                case 'Multiplicative noise'
                                    switch SplineType
                                        case 'LL'
                                            pp_mu = interp1(knots, par(1:M), 'linear', 'pp');
                                            pp_sigma = interp1(knots, par(M + 1:end), 'linear', 'pp');
                                            figure(1);
                                            plot(mesh_fine, ppval(pp_mu, mesh_fine), '-k');
                                            hold on;
                                            plot(knots, par(1:M), '*r');
                                            hold on;
                                            plot(knots, 0 .* knots, '-k');
                                            xlabel('State, x');
                                            ylabel('Drift function, \mu (x)');
                                            xlim([-inf inf]);
                                            ylim([-inf inf]);
                                            pos1 = get(gcf, 'Position');
                                            set(gcf, 'Position', pos1 - [pos1(3) / 2, 0, 0, 0]);
                                            figure(2);
                                            plot(mesh_fine, ppval(pp_sigma, mesh_fine), '-k');
                                            hold on;
                                            plot(knots, par(M + 1:end), '*r');
                                            xlabel('State, x');
                                            ylabel('Diffusion function, \sigma (x)');
                                            xlim([-inf inf]);
                                            ylim([0 inf]);
                                            pos2 = get(gcf, 'Position');
                                            set(gcf, 'Position', pos2 + [pos1(3) / 2, 0, 0, 0]);
                                            figure(3);
                                            A1 = ppval(pp_mu, mesh_fine);
                                            A2 = ppval(pp_sigma, mesh_fine);
                                            U = -2 .* cumtrapz(mesh_fine, A1 ./ A2.^2) + 2 .* log(A2);
                                            U = U - min(U) + 0.5;
                                            plot(mesh_fine, U, '-k');
                                            xlabel('State, x');
                                            ylabel('Effective potential, U_{eff}(x)');
                                            xlim([-inf inf]);
                                            ylim([0 inf]);
                                            pos3 = get(gcf, 'Position');
                                            set(gcf, 'Position', pos3 + [0, -504, 4, 0]);
                                        case {'CP', 'PP', 'CC', 'SCS', 'Approximate'}
                                            switch SplineType
                                                case 'CP'
                                                    pp_sigma = pchip(knots, par(M + 1:end));
                                                    pp_mu = spline(knots, par(1:M));
                                                case 'PP'
                                                    pp_sigma = pchip(knots, par(M + 1:end));
                                                    pp_mu = pchip(knots, par(1:M));
                                                case 'CC'
                                                    pp_sigma = spline(knots, par(M + 1:end));
                                                    pp_mu = spline(knots, par(1:M));
                                                case {'SCS', 'Approximate'}
                                                    pp_sigma = csaps(knots, par(M + 1:end));
                                                    pp_mu = csaps(knots, par(1:M));
                                            end
                                            figure(1);
                                            plot(mesh_fine, ppval(pp_mu, mesh_fine), '-k');
                                            hold on;
                                            plot(mesh_fine, 0 .* mesh_fine, '-k');
                                            hold on;
                                            plot(knots, par(1:M), '*r');
                                            xlabel('State, x');
                                            ylabel('Drift function, \mu (x)');
                                            xlim([-inf inf]);
                                            ylim([-inf inf]);
                                            pos1 = get(gcf, 'Position');
                                            set(gcf, 'Position', pos1 - [pos1(3) / 2, 0, 0, 0]);
                                            figure(2);
                                            plot(mesh_fine, ppval(pp_sigma, mesh_fine), '-k');
                                            hold on;
                                            plot(knots, par(M + 1:end), '*r');
                                            xlabel('State, x');
                                            ylabel('Diffusion function, \sigma (x)');
                                            xlim([-inf inf]);
                                            ylim([0 inf]);
                                            pos2 = get(gcf, 'Position');
                                            set(gcf, 'Position', pos2 + [pos1(3) / 2, 0, 0, 0]);
                                            figure(3);
                                            A1 = ppval(pp_mu, mesh_fine);
                                            A2 = ppval(pp_sigma, mesh_fine);
                                            U = -2 .* cumtrapz(mesh_fine, A1 ./ A2.^2) + 2 .* log(A2);
                                            U = U - min(U) + 0.5;
                                            plot(mesh_fine, U, '-k');
                                            xlabel('State, x');
                                            ylabel('Effective potential, U_{eff}(x)');
                                            xlim([-inf inf]);
                                            ylim([0 inf]);
                                            pos3 = get(gcf, 'Position');
                                            set(gcf, 'Position', pos3 + [0, -504, 4, 0]);
                                        case 'QQ'
                                            %                                         sp_mu=spapi(3,knots,par(1:M));sp_sigma=spapi(3,knots,par(M+1:end));
                                            pp_mu = spline(knots, par(1:M));
                                            pp_sigma = spline(knots, par(M + 1:end));
                                            figure(1);
                                            % plot(mesh_fine,fnval(sp_mu,mesh_fine),'-k');hold on;plot(mesh_fine,0 .*mesh_fine,'-k');hold on;plot(knots,par(1:M),'*r');
                                            plot(mesh_fine, ppval(pp_mu, mesh_fine), '-k');
                                            hold on;
                                            plot(mesh_fine, 0 .* mesh_fine, '-k');
                                            hold on;
                                            plot(knots, par(1:M), '*r');
                                            pos1 = get(gcf, 'Position');
                                            set(gcf, 'Position', pos1 - [pos1(3) / 2, 0, 0, 0]);
                                            xlabel('State, x');
                                            ylabel('Drift function, \mu (x)');
                                            xlim([-inf inf]);
                                            ylim([-inf inf]);
                                            figure(2);
                                            %                                         plot(mesh_fine,fnval(sp_sigma,mesh_fine),'-k');hold on;plot(knots,par(M+1:end),'*r');
                                            plot(mesh_fine, ppval(pp_sigma, mesh_fine), '-k');
                                            hold on;
                                            plot(knots, par(M + 1:end), '*r');
                                            xlabel('State, x');
                                            ylabel('Diffusion function, \sigma (x)');
                                            xlim([-inf inf]);
                                            ylim([0 inf]);
                                            pos2 = get(gcf, 'Position');
                                            set(gcf, 'Position', pos2 + [pos1(3) / 2, 0, 0, 0]);
                                            figure(3);
                                            %                                         A1=fnval(sp_mu,mesh_fine);A2=fnval(sp_sigma,mesh_fine);
                                            A1 = ppval(pp_mu, mesh_fine);
                                            A2 = ppval(pp_sigma, mesh_fine);
                                            U = -2 .* cumtrapz(mesh_fine, A1 ./ A2.^2) + 2 .* log(A2);
                                            U = U - min(U) + 0.5;
                                            plot(mesh_fine, U, '-k');
                                            xlabel('State, x');
                                            ylabel('Effective potential, U_{eff}(x)');
                                            xlim([-inf inf]);
                                            ylim([0 inf]);
                                            pos3 = get(gcf, 'Position');
                                            set(gcf, 'Position', pos3 + [0, -505, 4, 0]);
                                    end
                            end
                            drawnow;
                        end
                    end
                end

                if fitness > Alpha_score && fitness < Beta_score
                    Beta_score = fitness; % Update beta
                    Beta_pos = Positions(i, :);
                end

                if fitness > Alpha_score && fitness > Beta_score && fitness < Delta_score
                    Delta_score = fitness; % Update delta
                    Delta_pos = Positions(i, :);
                end
            end

            %% Calculate the candiadate position Xi-GWO
            a = 2 - iter * ((2) / Max_iter); % a decreases linearly from 2 to 0

            % Update the Position of search agents including omegas


            %     for i=1:N
            %         for j=1:dim
            %             r1=rand(); % r1 is a random number in [0,1]
            %             r2=rand(); % r2 is a random number in [0,1]
            %
            %             A1=2*a*r1-a;                                    % Equation (3.3)
            %             C1=2*r2;                                        % Equation (3.4)
            %
            %             D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j));    % Equation (3.5)-part 1
            %             X1=Alpha_pos(j)-A1*D_alpha;                     % Equation (3.6)-part 1
            %
            %             r1=rand();
            %             r2=rand();
            %
            %             A2=2*a*r1-a;                                    % Equation (3.3)
            %             C2=2*r2;                                        % Equation (3.4)
            %
            %             D_beta=abs(C2*Beta_pos(j)-Positions(i,j));      % Equation (3.5)-part 2
            %             X2=Beta_pos(j)-A2*D_beta;                       % Equation (3.6)-part 2
            %
            %             r1=rand();
            %             r2=rand();
            %
            %             A3=2*a*r1-a;                                    % Equation (3.3)
            %             C3=2*r2;                                        % Equation (3.4)
            %
            %             D_delta=abs(C3*Delta_pos(j)-Positions(i,j));    % Equation (3.5)-part 3
            %             X3=Delta_pos(j)-A3*D_delta;                     % Equation (3.5)-part 3
            %
            %             X_GWO(i,j)=(X1+X2+X3)/3;                        % Equation (3.7)
            %         end
            %         X_GWO(i,:) = boundConstraint(X_GWO(i,:), Positions(i,:), lu);
            %         Fit_GWO(i) = fobj(X_GWO(i,:));
            %     end

            %**********My Additions****************************
            r1 = rand(N, dim);
            r2 = rand(N, dim);
            A1 = 2 .* a .* r1 - a;
            C1 = 2 .* r2;
            Alphapos = repmat(Alpha_pos, N, 1);
            D_alpha = abs(C1 .* Alphapos - Positions);
            X1 = Alphapos - A1 .* D_alpha;

            r1 = rand(N, dim);
            r2 = rand(N, dim);
            A2 = 2 .* a .* r1 - a;
            C2 = 2 .* r2;
            Betapos = repmat(Beta_pos, N, 1);
            D_beta = abs(C2 .* Betapos - Positions);
            X2 = Betapos - A2 .* D_beta;

            r1 = rand(N, dim);
            r2 = rand(N, dim);
            A3 = 2 .* a .* r1 - a;
            C3 = 2 .* r2;
            Deltapos = repmat(Delta_pos, N, 1);
            D_delta = abs(C3 .* Deltapos - Positions);
            X3 = Deltapos - A3 .* D_delta;

            X_GWO = (X1 + X2 + X3) / 3;
            IDX = X_GWO(:, M + 1:end) <= 0;
            IDX = [false(N, M) IDX];
            X_GWO(IDX) = (Positions(IDX) + LB(IDX)) ./ 2; % We only correct for the lowe bounds of diffusion function when they get negative or 0

            parfor (i = 1:N, numWorkers)
                fobj1 = fobj;
                Fit_GWO(i) = fobj1(X_GWO(i, :));
            end
            %***********************************************************************************
            %% Calculate the candiadate position Xi-DLH
            radius = pdist2(Positions, X_GWO, 'euclidean'); % Equation (10)
            dist_Position = squareform(pdist(Positions));
            r1 = randperm(N, N);

            X_DLH = zeros(N, dim);
            for t = 1:N
                neighbor(t, :) = (dist_Position(t, :) <= radius(t, t));
                [~, Idx] = find(neighbor(t, :) == 1); % Equation (11)
                random_Idx_neighbor = randi(size(Idx, 2), 1, dim);

                %         for d=1:dim
                %             X_DLH(t,d) = Positions(t,d) + rand .*(Positions(Idx(random_Idx_neighbor(d)),d)- Positions(r1(t),d));                      % Equation (12)
                %         end
                R = rand(1, dim);
                X_DLH(t, :) = Positions(t, :) + R .* (Positions(sub2ind(size(Positions), Idx(random_Idx_neighbor), 1:dim)) - Positions(r1(t), :)); %My addition

                %         X_DLH(t,:) = boundConstraint(X_DLH(t,:), Positions(t,:), lu;

                idx = X_DLH(t, M + 1:end) <= 0;
                idx = [false(1, M) idx]; %My addition
                X_DLH(t, idx) = (Positions(t, idx) + lb(idx)) ./ 2; %My addition

                %         Fit_DLH(t) = fobj(X_DLH(t,:));
            end

            parfor (t = 1:N, numWorkers) %My addition
                fobj1 = fobj;
                Fit_DLH(t) = fobj1(X_DLH(t, :));
            end

            %% Selection
            tmp = Fit_GWO < Fit_DLH; % Equation (13)
            tmp_rep = repmat(tmp', 1, dim);

            tmpFit = tmp .* Fit_GWO + (1 - tmp) .* Fit_DLH;
            tmpPositions = tmp_rep .* X_GWO + (1 - tmp_rep) .* X_DLH;

            %% Updating
            tmp = pBestScore <= tmpFit; % Equation (13)
            tmp_rep = repmat(tmp', 1, dim);

            pBestScore = tmp .* pBestScore + (1 - tmp) .* tmpFit;
            pBest = tmp_rep .* pBest + (1 - tmp_rep) .* tmpPositions;

            Fit = pBestScore;
            Positions = pBest;

            %%
            iter = iter + 1;
            neighbor = zeros(N, N);
            %     Convergence_curve(iter) = Alpha_score;
        end
    end
    if ~isempty(outputFcn)
        outputFcn(optimValues, 'done');
    end
end




