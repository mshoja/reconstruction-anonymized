
function [lb, ub, Leg_Points] = LegitimateRegion(cost, dim, lb, ub, index_sigma)
%     if useparallel
%         numworkers=Inf;
%     else
%         numworkers=1;
%     end

Npoints = 30; %By the Central Limit Theorem N = 30 points should be good (for the purpose of finding the center of mass
% of legitimate points). A very great acuracy is not needed here

Leg_Points = GenerateLegitimatePoints1(Npoints, dim, lb, ub, index_sigma, cost);

%Estimating the center of mass for the obtained legitimate points
fprintf('\nEstimating the center of mass for the already found legitimate solutions\n');

d = abs(Leg_Points(:, dim + 1) - Leg_Points(1, dim + 1));
d_max = max(d);
w = -d ./ d_max + 1;
w = w(:);
w = w ./ sum(w); %weights
CM = sum(Leg_Points(:, 1:dim) .* ones(1, dim) .* w); %Center of mass

CM


% Roughly finging the radius of a box, centered at CM, under which 25 percent of points on the surface of this box are legitimate
fprintf('\nEstimating the radius of a box centered at the center of mass of legitimate points which contains 25%% legitimate points on its surface\n');
Npoints = 10^2; % By the Central Limit Theorem 100 points should be great (we chose more than 30 since here it is not really expensive)
p = 0.25; % p is the percentage of legitimate points we wish to exist on the surface of the considered box
cost_radius = @(r)Cost_radius(r,p,CM,Npoints,cost);
r0 = max(abs(CM));
r25 = r0;
lb1 = 0;
ub1 = 10^2 * r0; % The upper limit of THIS optimization can be as big as possible in theory. However, since the objective function is extremely flat for very big r we are limitted to
% a smaller upper bound. Fortunately, we really do not need to go way beyond r0
% options = optimoptions('patternsearch', 'Display', 'none');
options = optimoptions('fmincon', 'Display', 'none');


% 25% radius
while abs(r25 - r0) < 10^(-3) || abs(ub1 - r25) < 10^(-3) % The first condition is violated if the chosen upper limit, ub1, is extremely big. Since the objective function becomes very flat for big radius r then
    % in this case r25=r0 is the outcome. The second condition is violated when ub1 is very small in which case the optimization algorithm cannot go beyond ub1 and hit it so that the  outcome of optimization is r=ub1 (almost often the second condition will not be violated. We never experienced this)
    % [r25, ~] = patternsearch(cost_radius, r0, [], [], [], [], lb1, ub1, [], options);
    [r25, ~] = fmincon(cost_radius,r0,[],[],[],[],lb1,ub1,[],options);
    if abs(r25 - r0) < 10^(-1) %this is a sign for the failure of THIS optimization since upper limit is so big
        ub1 = ub1 / 2;
    end
    if abs(ub1 - r25) < 10^(-3) %this is a sign for the failure of THIS optimization since upper limit is so big
        ub1 = 2 * ub1;
    end
end
lb = CM - r25;
ub = CM + r25;
end

function c = Cost_radius(r, p, CM, Npoints, cost)
lb = CM - r;
ub = CM + r;
RS = RandomOnBoxSurface(lb, ub, Npoints);

A = zeros(1, Npoints);
parfor i = 1:Npoints
    cost1 = cost;
    A(i) = cost1(RS(i, :));
end
A(A == realmax) = nan;
ratio = nnz(~isnan(A)) / Npoints;
c = (ratio - p)^2;
% c=abs(ratio-p);
end

function RS = RandomOnBoxSurface(lb, ub, N)

% This code generates N random points on the surface of the box being
% defined by a vetor of lower bounds lb and a vector of upper bounds ub

dim = length(lb);

lb = lb(:).';
ub = ub(:).';
LB = repmat(lb, [N 1]);
UB = repmat(ub, [N 1]);

%Generate N random points inside the unit box of dimension (dim-1). We follow sobol sampling
p1 = sobolset(dim - 1, 'Skip', 1e3, 'Leap', 1e2);
p1 = scramble(p1, 'MatousekAffineOwen');
RS1 = net(p1, N);
% RS1=rand(N,dim-1);   % if you wish to follow simple random sampling!

I = rand(N, 1) > 0.5; %Generate a vector of N binary elements (each 0 or 1 have 50% chance)
RS = [RS1 I]; %place the binary vector at the last column (or any other column)

% Shuffle the rows of A
[~, idx] = sort(rand(N, dim), 2);
idx = (idx - 1) * N + ndgrid(1:N, 1:dim);
RS = RS(idx); %each row of B is a point on the surface of unit box of dimension dim

%Translate points (on the surface of unit box) to the surface of the box defined by the lower vector lb and upper vector ub
RS = LB + (UB - LB) .* RS;
end

function Leg_Points = GenerateLegitimatePoints1(N, dim, lb, ub, ndx_sigma, cost)

% A short explanation on what this code does
% Because it is not clear 1) how big the searching space should be 2) if we choose very smal searching region we find 
% no legitimate solution 3) if we choose rather big searching space we find legitimate solutions but it might be computationaly
% expensive,  we therefore consider 4 searching regions and randomly search within them.

% Note : the above mentioned 4 searching regions are co-centric. However, when the problem dimensions gets bigger 

ran = (ub(1)-lb(1))./2;
x = ub-ran;  % this is the Euler solution

lb1 = lb;ub1 = ub;
lb2 = x-ran/2;ub2 = x+ran/2;
lb3 = x-ran/4;ub3 = x+ran/4;
lb4 = x-ran/8;ub4 = x+ran/8;

ndx_sigma1 = ndx_sigma & lb1 < eps;
ndx_sigma2 = ndx_sigma & lb2 < eps;
ndx_sigma3 = ndx_sigma & lb3 < eps;
ndx_sigma4 = ndx_sigma & lb4 < eps;

lb1(ndx_sigma1) = eps;
lb2(ndx_sigma2) = eps;
lb3(ndx_sigma3) = eps;
lb4(ndx_sigma4) = eps;  

Leg_Points = search(lb1, ub1, lb2, ub2, lb3, ub3, lb4, ub4, N, dim, cost);


for i=1:10
    i
    a=search(lb1, ub1, lb2, ub2, lb3, ub3, lb4, ub4, N, dim, cost);

end


fprintf('A first legitimate solution is already found\n');

% We find a uni-dimensinal line segment along which we can move forward and backward starting from each
% legitimate points. This helps a lot to know the area (= box) within which majority of legitimate points live
% Leg_Points = Leg_Points(:,1:dim);
segment_forward = zeros(1,dim);
segment_backward = zeros(1,dim);
ub1=nan(1,dim);
lb1=nan(1,dim);
tol = 10^(-2);
for i=1:size(Leg_Points,1)
    point = Leg_Points(i,1:dim);
    parfor j=1:dim
        cost1 = @(r)Cost_segment(r, point, j , cost);
        segment_forward(j) = bisection_forward(cost1, 0, ub(j)-point(j), tol);
        segment_backward(j) = bisection_backward(cost1, lb(j)-point(j), 0, tol);
    end
    ub1 = max(ub1, point+segment_forward);
    lb1 = min(lb1, point+segment_backward);
end

% In below we try to enlarge the box [lb1, ub1] through branching (branching = we fix all dimensions but only vary one
% of them in forward and backward directions from a legitimate solution as far as this is possible). This is an 
% iterative process and the more one does this the better (we consider 10 iterations).

i = 1;
while i<11
    point = lb1+(ub1-lb1).*rand(1, dim);
    if cost(point)<realmax
        parfor j=1:dim
            cost1 = @(r)Cost_segment(r, point, j , cost);
            segment_forward(j) = bisection_forward(cost1, 0, ub(j)-point(j), tol);
            segment_backward(j) = bisection_backward(cost1, lb(j)-point(j), 0, tol);
        end
        ub1 = max(ub1, point+segment_forward);
        lb1 = min(lb1, point+segment_backward);
        i =i+1;
    end
end

% % In bellow we try to find N legitimate points within the box [lb1, ub1]
tic;
while size(Leg_Points, 1)<N
    size(Leg_Points, 1)
    parfor i=1:10
        cost1 = cost;
        point = lb1+(ub1-lb1).*rand(1, dim);
        f = cost1(point);
        if f<realmax
            Leg_Points = [Leg_Points; [point f]];
        end
    end  
end
toc;

% tic;
% while size(Leg_Points, 1)<N
%     point = lb1+(ub1-lb1).*rand(1, dim);
%     f = cost(point);
%     if f<realmax
%         Leg_Points = [Leg_Points; [point f]];
%     end   
% end
% toc;
Leg_Points = sortrows(Leg_Points, dim + 1);
disp([num2str(size(Leg_Points, 1)) ' legitimate solutions are already found']);
fprintf('\n');

% vpa(Leg_Points,8)
end

function c = Cost_segment(r, point, axis , cost)
point(axis) = point(axis)+r;
f = cost(point);
c = 1;
if f==realmax
    c = 0;
end
end

function x3 =bisection_forward(f, lb, ub, tol)
x1 = lb;
x2 = ub;
if f(x2)
    x3 = x2;
else
    d = 2*tol;
    while d>tol
        x3 = (x1+x2)/2;
        f3 = f(x3);
        if f3
            x1=x3;
        else
            x2 = x3;
        end
        d = x2-x1;
    end
end
idx = exist('f3','var');
if idx && ~f3
    x3 = x1;
end
end

function x3 =bisection_backward(f, lb, ub, tol)
x1 = lb;
x2 = ub;
if f(x1)
    x3 = x1;
else
    d = 2*tol;
    while d>tol
        x3 = (x1+x2)/2;
        f3 = f(x3);
        if f3
            x2=x3;
        else
            x1 = x3;
        end
        d = x2-x1;
    end
end
idx = exist('f3','var');
if idx && ~f3
    x3 = x2;
end
end

function Leg_Points = search(lb1, ub1, lb2, ub2, lb3, ub3, lb4, ub4, N, dim, cost)
Leg_Points = [zeros(N, dim) realmax .* ones(N, 1)];
n1 = 5 * N;
LB1 = repmat(lb1, [n1 1]);UB1 = repmat(ub1, [n1 1]);
LB2 = repmat(lb2, [n1 1]);UB2 = repmat(ub2, [n1 1]);
LB3 = repmat(lb3, [n1 1]);UB3 = repmat(ub3, [n1 1]);
LB4 = repmat(lb4, [n1 1]);UB4 = repmat(ub4, [n1 1]);

n = 0;
while n == 0
    m = randi(4);
    switch m
        case 1
            LB = LB1;UB = UB1;
        case 2
            LB = LB2;UB = UB2;
        case 3
            LB = LB3;UB = UB3;
        case 4
            LB = LB4;UB = UB4;
    end

    % Sampling based on sobol set (this is better than random and stratified random samplings especially in higher dimensions )
    p = sobolset(dim, 'Skip', 1e3, 'Leap', 1e2);
    p = scramble(p, 'MatousekAffineOwen');
    RS = net(p, n1); % This is a random sample whose entities are in (0,1)
    RS = LB + (UB - LB) .* RS;
    f = zeros(n1, 1);
    parfor i = 1:n1
        cost1 = cost;
        f(i) = cost(RS(i, :));
    end
    Leg_Points = [RS f; Leg_Points];
    Leg_Points = sortrows(Leg_Points, dim + 1);
    b = Leg_Points(:, end);
    idx = b ~= realmax;
    Leg_Points = Leg_Points(idx, :);
    n = size(Leg_Points,1);
end
end



%In case we do not like to use patternsearch optimization to find the
%radius of searching region we need the following code
%************************************************************************************************************************************************
%**********************************JUST IN CASE WE DO NOT LIKE TO USE PATTERNSEARCH**************************************************************
%************************************************************************************************************************************************

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
function [Alpha_score, Alpha_pos] = IGWO(dim, N, Max_iter, lb, ub, fobj)


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
% Positions = boundConstraint (Positions, Positions, lu);

%My addition
%Simple random sampling

RS = rand(N, dim);
LB = repmat(lb, [N 1]);
UB = repmat(ub, [N 1]);
Positions = LB + (UB - LB) .* RS;

% Calculate objective function for each wolf
parfor i = 1:size(Positions, 1)
    fobj1 = fobj;
    Fit(i) = fobj1(Positions(i, :));
end

% Personal best fitness and position obtained by each wolf
pBestScore = Fit;
pBest = Positions;

neighbor = zeros(N, N);
% Convergence_curve=zeros(1,Max_iter);
iter = 0; % Loop counter

% Main loop
Objective = realmax;
% while iter < Max_iter
while Objective > 10^(-3)
    for i = 1:size(Positions, 1)
        fitness = Fit(i);

        % Update Alpha, Beta, and Delta
        if fitness < Alpha_score
            Alpha_score = fitness; % Update alpha
            Alpha_pos = Positions(i, :);
            Objective = Alpha_score;
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

    % Calculate the candiadate position Xi-GWO
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
    %
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
    idx_LB = X_GWO < LB;
    idx_UB = X_GWO > UB;

    % Correcting for the violation, if any, of bound constraints
    X_GWO(idx_LB) = (Positions(idx_LB) + LB(idx_LB)) ./ 2;
    X_GWO(idx_UB) = (Positions(idx_UB) + UB(idx_UB)) ./ 2;

    parfor i = 1:N
        fobj1 = fobj;
        Fit_GWO(i) = fobj1(X_GWO(i));
    end

    % Calculate the candiadate position Xi-DLH
    radius = pdist2(Positions, X_GWO, 'euclidean'); % Equation (10)
    dist_Position = squareform(pdist(Positions));
    r1 = randperm(N, N);

    for t = 1:N
        neighbor(t, :) = (dist_Position(t, :) <= radius(t, t));
        [~, Idx] = find(neighbor(t, :) == 1); % Equation (11)
        random_Idx_neighbor = randi(size(Idx, 2), 1, dim);

        %         for d=1:dim
        %             X_DLH(t,d) = Positions(t,d) + rand .*(Positions(Idx(random_Idx_neighbor(d)),d)...
        %                 - Positions(r1(t),d));                      % Equation (12)
        %         end
        %         X_DLH(t,:) = boundConstraint(X_DLH(t,:), Positions(t,:), lu);

        R = rand(1, dim);
        X_DLH(t, :) = Positions(t, :) + R .* (Positions(sub2ind(size(Positions), Idx(random_Idx_neighbor), 1:dim)) - Positions(r1(t), :)); %My addition
        idx_lb = X_DLH(t, :) < lb;
        idx_ub = X_DLH(t, :) > ub; %My addition
        X_DLH(t, idx_lb) = (Positions(t, idx_lb) + lb(idx_lb)) ./ 2; %My addition
        X_DLH(t, idx_ub) = (Positions(t, idx_ub) + ub(idx_ub)) ./ 2; %My addition

        Fit_DLH(t) = fobj(X_DLH(t, :));
    end

    % Selection
    tmp = Fit_GWO < Fit_DLH; % Equation (13)
    tmp_rep = repmat(tmp', 1, dim);

    tmpFit = tmp .* Fit_GWO + (1 - tmp) .* Fit_DLH;
    tmpPositions = tmp_rep .* X_GWO + (1 - tmp_rep) .* X_DLH;

    % Updating
    tmp = pBestScore <= tmpFit; % Equation (13)
    tmp_rep = repmat(tmp', 1, dim);

    pBestScore = tmp .* pBestScore + (1 - tmp) .* tmpFit;
    pBest = tmp_rep .* pBest + (1 - tmp_rep) .* tmpPositions;

    Fit = pBestScore;
    Positions = pBest;

    %
    iter = iter + 1;
    neighbor = zeros(N, N);
    %     Convergence_curve(iter) = Alpha_score;
end
end














