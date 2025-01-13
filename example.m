
%%Example data set Ornstad Ulenberg model
%we load the data OUdata1D and define




%load the data
S=load('OUdata1D');
data1=S.data(1:10000);
%This data set is created the following mu and sigma functions 
mu = @(x,par)par(1)*x;
sigma = @(x,par)par(2)+zeros(size(x));
%the real parameters of the model
realpar = [-1 1];
%the mu and sigma on which the data is based (used for plotting later)
xtrue = linspace(min(data1), max(data1), 200);
mutrue = mu(xtrue, realpar);
sigmatrue = sigma(xtrue, realpar);

%Relaxation time, fit to autocorrelation function: acf=exp(-c*lag)
%Relaxation time is 1/c
fprintf('The relaxation time = %g\n', RelaxationTime(data1))

%%We assume additive noise and use cubic splines to describe the drift
%%function (mu)
%optionally you can use a number of different starts

% result1 = euler_reconstruction(data1, dt, 'nKnots', [4 1], 'spline', 'CC',   ...
%     'lb', [zeros(1, 4) - 10, zeros(1, 1)], 'ub', zeros(1, 4 + 1) + 10,  ...
%     'solver', 'gwo', 'search_agents', 5);

% result1 = euler_reconstruction(data1, dt, 'nKnots', [2 2], 'spline', 'CC',  ...
%     'lb', [zeros(1, 2) - 10, zeros(1, 2)+eps], 'ub', zeros(1, 2 + 2) + 10,  ...
%     'solver', 'gwo', 'search_agents', 5);

L=-2;R=1.8;
result1 = euler_reconstruction(data1, dt, 'nKnots', [6 6], 'spline', 'CC', 'L', L, 'R', R,  ...
    'lb', [zeros(1, 6) - 10, zeros(1, 6)], 'ub', zeros(1, 6 + 6) + 10,  ...
    'solver', 'gwo', 'search_agents', 5);
%% 

%plot the Euler results
% plot_results(result1, xtrue, mutrue, sigmatrue)

%find parameter range with many legitimate points >25%
% this step is only needed for the fmincon solver
% we use the results of the Euler +/- the prev_range as initial lower and
% upper bounds

legpoints = legitimate_points(data1, dt, 'prev', result1, 'prev_range', 1, 'j', 3, 'k', 5);

save legpoints legpoints;
%save these results for future use

load legpoints

%run the hermite reconstruction using the fmincon solver, with 5 random initial
%conditions
%
result2 = hermite_reconstruction(data1, dt, 'prev', legpoints, 'solver', 'fmincon', 'search_agents', 5, 'UseParallel', true);
%plot the hermite results
plot_results(result2, xtrue, mutrue, sigmatrue);


%% parametric reconstruction 
%if we know or assume a model we can also fit the parameters of both mu and
%sigma. Instead of a spline we need to supply a function for mu and sigma
%in the form (where par is a vector of all parameters both of mu and sigma.
%for example:
%mu=@(x,par)par(1).*x.^2+par(2).*x

%or
mu = @(x,par)par(1)*x;
sigma = @(x,par)par(2)+zeros(size(x));
%First the euler method:

result3 = euler_reconstruction(data1, dt, 'mu', mu, 'sigma', sigma, 'lb', [-200 eps], 'ub', [200 200]);

%we plot the results
plot_results(result3, xtrue, mutrue, sigmatrue);

%again if we wat to use Hermite reconstruction with fmincon, we first need
%to find legitimate points:
legpoints2 = legitimate_points(data1, dt, 'prev', result3, 'prev_range', 0.5, 'j', 3, 'k', 5);

%using the legitimate points we continue:
result4 = hermite_reconstruction(data1, dt, 'prev', legpoints2, 'solver', 'gwo', 'search_agents', 10);

%plot the results
plot_results(result4, xtrue, mutrue, sigmatrue);
