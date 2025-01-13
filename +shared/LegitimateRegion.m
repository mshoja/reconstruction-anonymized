
function Leg_Points = LegitimateRegion(cost, dim, lb, ub, search_agents, index_sigma)

Npoints = search_agents;
Leg_Points = GenerateLegitimatePoints(Npoints, dim, lb, ub, index_sigma, cost);

Leg_Points = sortrows(Leg_Points, dim+1);
end

% A short explanation on what this code does
% Because it is not clear 1) how big the searching space should be 2) if we choose very smal searching region we find 
% no legitimate solution 3) if we choose rather big searching space we find legitimate solutions but it might be computationaly
% expensive,  we therefore consider the following searching algorithm which updates the searching region based on 
% its previous attemts (i.e., if previous tries were successful it follows them, otherwise it will not. But, to make 
% sure algorithm always works there is also a random component to this strategy).

function Leg_Points = GenerateLegitimatePoints(N, dim, lb, ub, ndx_sigma, cost)

% Leg_Points = [zeros(N, dim) realmax .* ones(N, 1)];
% LB = repmat(lb, [N 1]);UB = repmat(ub, [N 1]);n = 0;
% while n<N
%     p = sobolset(dim, 'Skip', 1e3, 'Leap', 1e2);
%     p = scramble(p, 'MatousekAffineOwen');
%     RS = net(p, N); % This is a random sample whose entities are in (0,1)
%     RS = LB + (UB - LB) .* RS;
%     f = zeros(N, 1);
%     parfor i = 1:N
%         cost1 = cost;
%         f(i) = cost1(RS(i, :));
%     end
%     Leg_Points = [RS f; Leg_Points];
%     Leg_Points = sortrows(Leg_Points, dim + 1);
%     b = Leg_Points(:, end);
%     idx = b ~= realmax;
%     Leg_Points = Leg_Points(idx, :);
%     n = size(Leg_Points,1);
% end


ran = (max(ub) - min(lb))/2;
ndx_sigma = ndx_sigma & lb < eps;
lb(ndx_sigma) = eps;
Leg_Points = [];
if isfile('population.csv')
    population = readmatrix('population.csv');
else
    population = [];
end
for i=1:N
    [LP, population] = search(population, ran, ub, ndx_sigma, N, dim, cost);
    Leg_Points = [Leg_Points; LP];
    if i==1
        fprintf('\nA first legitimate solution is already found\n');
    else
        disp([num2str(i) ' legitimate solutions are already found']);
    end
end
disp([num2str(size(Leg_Points, 1)) ' legitimate solutions are found']);
writematrix(unique(population),'population.csv')
end

function [Leg_Points, population] = search(population, ran, ub, ndx_sigma, N, dim, cost)
Leg_Points = [zeros(N, dim) realmax .* ones(N, 1)];
% n1=10;
n = 0;
while n == 0
%     r0 = ran*rand;
    a = rand;
    r0 = ran*a;
    r = randsample([population r0],1);
    X = ub-r;   % if rand = 1 then X = ub-ran = Euler and we search in the box centered at Euler with the biggest side
    lb = X-r;
    ub = X+r;
    ndx_sigma = ndx_sigma & lb < eps;
    lb(ndx_sigma) = eps;

    n1 = 10+ceil(a/0.01);
    LB = repmat(lb, [n1 1]);UB = repmat(ub, [n1 1]);

    % Sampling based on sobol set (this is better than random and stratified random samplings especially in higher dimensions )
    p = sobolset(dim, 'Skip', 1e3, 'Leap', 1e2);
    p = scramble(p, 'MatousekAffineOwen');
    RS = net(p, n1); % This is a random sample whose entities are in (0,1)
    RS = LB + (UB - LB) .* RS;
    f = zeros(n1, 1);
    parfor i = 1:n1
        cost1 = cost;
        f(i) = cost1(RS(i, :));
    end
    Leg_Points = [RS f; Leg_Points];
    Leg_Points = sortrows(Leg_Points, dim + 1);
    b = Leg_Points(:, end);
    idx = b ~= realmax;
    Leg_Points = Leg_Points(idx, :);
    n = size(Leg_Points,1);
end
population = [population r];
end