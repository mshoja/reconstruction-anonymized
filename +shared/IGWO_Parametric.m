%___________________________________________________________________%
%  An Improved Grey Wolf Optimizer for Solving Engineering          %
%  Problems (I-GWO) source codes version 1.0                        %
%                                                                   %
%  Developed in MATLAB R2018a                                       %
%                                                                   %
%  Author and programmer: M. H. Nadimi-Shahraki, S. Taghian, S. Mirjalili                                           %
%                                                                   %
% e-Mail: nadimi@ieee.org, shokooh.taghian94@gmail.com,                                         ali.mirjalili@gmail.com                                   %
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
function [Alpha_pos, Alpha_score] = IGWO_Parametric(dim, N, Max_iter, lb, ub, fobj, M, UseParallel, outputFcn)
    if nargin < 9
        outputFcn = [];
    end
    if UseParallel
        numWorkers = Inf;
    else
        numWorkers = 1;
    end %My addition


    % lu = [lb .* ones(1, dim); ub .* ones(1, dim)];

    % Initialize alpha, beta, and delta positions
    Alpha_pos = zeros(1, dim);
    Alpha_score = inf; %change this to -inf for maximization problems

    Beta_pos = zeros(1, dim);
    Beta_score = inf; %change this to -inf for maximization problems

    Delta_pos = zeros(1, dim);
    Delta_score = inf; %change this to -inf for maximization problems

    % Initialize the positions of wolves
    % Positions=initialization(N,dim,ub,lb);
    %**********************************
    %My additions
    optimValues = struct('bestfval', NaN, 'bestx', zeros(size(lb)), 'iteration', 0);
    if ~isempty(outputFcn)
        outputFcn(optimValues, 'init');
    end
    LB = repmat(lb, [N 1]);
    UB = repmat(ub, [N 1]);
    Positions = LB + (UB - LB) .* rand(N, dim);
    %*********************************
    % Positions = boundConstraint (Positions, Positions, lu);

    % Calculate objective function for each wolf
    for i = 1:size(Positions, 1)
        Fit(i) = fobj(Positions(i, :));
    end

    % Personal best fitness and position obtained by each wolf
    pBestScore = Fit;
    pBest = Positions;

    neighbor = zeros(N, N);
    % Convergence_curve=zeros(1,Max_iter);
    iter = 0; % Loop counter

    %% Main loop
    stop=false;
    while iter < Max_iter&&~stop
        for i = 1:size(Positions, 1)
            fitness = Fit(i);
     
            % Update Alpha, Beta, and Delta
            if fitness < Alpha_score
                Alpha_score = fitness; % Update alpha
                Alpha_pos = Positions(i, :);

                %My additions 
                par = Alpha_pos;
                if nargin == 7 && length(par(M:end)) == 1 %this is the case of additive Euler reconstruction 
                    par(M:end) = abs(par(M:end)); %In Euler reconstruction the objective function is invariant under negation of diffusion vector: any legitimate (i.e., positive) vector of diffusion function corresponds with an illegitimate (i.e., negative) solution whose absolute vale is the same
                end
                optimValues = struct('bestfval', Alpha_score, 'bestx', par, 'iteration', iter);
                if ~isempty(outputFcn)
                    stop = outputFcn(optimValues, 'iter');
                else
                    if isequal(Max_iter, realmax) %Displaying the results in the command window
                        if isequal(Alpha_score, realmax)
                            fprintf('\nThe algorithm could not find a legitimate solution so far. If this lasts long then you should stop the code, pay attention to the following possible reasons for the failure, and re-run again\n');
                            fprintf('\n1) There is a possibility that you have considered a very big (more probable) or very small (less probable) bound around the Euler solution. You should shrink (or expand) the bounds around the Euler solution');
                            fprintf('\n2) If your data does not have a high-resolution then you do not need to use a very small ''K'' value which might fail (even if it works is not useful for your data)');
                            fprintf('\nIf (1),(2) are not the reason for the failure then consider the following actions:');
                            fprintf('\n1) Consider a simpler model (e.g., of polynomil type)');
                            fprintf('\n2) A nice option is to consider spline models, instead\n\n');
                        else
                            disp(['Iteration = ', num2str( iter)]);
                            disp('Estimated parameters : ');
                            disp(num2str(par))
                            disp(['Approximate value of objective function (negative sum of log-likelihoods) : ' num2str(Alpha_score)]);
                        end
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
        %     for i=1:size(Positions,1)
        %         for j=1:size(Positions,2)
        %             
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
        % %         X_GWO(i,:) = boundConstraint(X_GWO(i,:), Positions(i,:), lu);
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
        IDX_LB = X_GWO < LB;
        IDX_UB = X_GWO > UB;
        X_GWO(IDX_LB) = (Positions(IDX_LB) + LB(IDX_LB)) ./ 2; %Correcting for lower bounds
        X_GWO(IDX_UB) = (Positions(IDX_UB) + LB(IDX_UB)) ./ 2; %Correcting for upper bounds


        parfor (i = 1:N, numWorkers)
            fobj1 = fobj;
            Fit_GWO(i) = fobj1(X_GWO(i, :));
        end
        %*************************************************************
        %% Calculate the candiadate position Xi-DLH
        radius = pdist2(Positions, X_GWO, 'euclidean'); % Equation (10)
        dist_Position = squareform(pdist(Positions));
        r1 = randperm(N, N);
     
        X_DLH = zeros(N, dim); %My addition
        for t = 1:N
            neighbor(t, :) = (dist_Position(t, :) <= radius(t, t));
            [~, Idx] = find(neighbor(t, :) == 1); % Equation (11)             
            random_Idx_neighbor = randi(size(Idx, 2), 1, dim);
     
            %         for d=1:dim
            %             X_DLH(t,d) = Positions(t,d) + rand .*(Positions(Idx(random_Idx_neighbor(d)),d)...
            %                 - Positions(r1(t),d));                      % Equation (12)
            %         end
            %         X_DLH(t,:) = boundConstraint(X_DLH(t,:), Positions(t,:), lu);

            %My additions
            R = rand(1, dim);
            X_DLH(t, :) = Positions(t, :) + R .* (Positions(sub2ind(size(Positions), Idx(random_Idx_neighbor), 1:dim)) - Positions(r1(t), :));
     
            idx_lb = X_DLH(t, :) < lb;
            idx_ub = X_DLH(t, :) > ub;
            X_DLH(t, idx_lb) = (Positions(t, idx_lb) + lb(idx_lb)) / 2; %Correcting for lower bounds 
            X_DLH(t, idx_ub) = (Positions(t, idx_ub) + ub(idx_ub)) / 2; %Correcting for upper bounds
        end
     
        %My addition
        parfor (t = 1:N, numWorkers)
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




