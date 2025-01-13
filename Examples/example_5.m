%In this example we analyse low resolution downsampled data set
%generated with the linear Ornstein-Uhlenbeck (OU) model
%Here we apply the spline method.

S = load('OUdata1D.mat');
data = S.data; %load the data
%We downsample to illustrate the performance for low resolution data
data = data(1:100:end); % Only every 100 data points are considered

R = RelaxationTime(data);

fprintf('Data resolution\nRelaxation time = %g\n', R);
if R > 100
    disp('High resolution data')
elseif R > 50
    disp('Medium resolution data')
elseif R > 1
    disp('Low resolution data')
else
    disp('Extremely low resolution data');
end
%When applying a  spline  model it is very important to consider shrinking the state space
%Typically, the range of state space is, by default, [min(data) max(data)] 
%Here= that is [-2.6289 2.9383].
%However, the dataset contains only 2 data points larger than 2.5 and 2 data points 
%smaller than -2.5 out of a total of 10,000 data points. Therefore it is strongly recommended to
%shrink the range of the state space 
L = -2.5;
R = 2.5;
%We assign ‘nan’ to those few data points falling outside this range.
data(data < L | data > R) = nan;
dt = 1;
nmu = 7; %number of knots for mu
nsigma = 7; %number of knots for sigma
%approximate Euler analysis
result11 = euler_reconstruction(data, dt, 'nKnots', [nmu nsigma], 'spline', 'QQ', 'L', L, 'R', R, ..., 
    'lb', [zeros(1, nmu) - 10, zeros(1, nsigma) + eps], 'ub', zeros(1, nmu + nsigma) + 10, 'solver', 'fmincon', 'search_agents', 5);
%find legitimate starting points based on the Euler solution
disp('searching legitimate points');
legpoints = legitimate_points(data, dt, 'prev', result11, 'prev_range', 0.5, 'j', 3, 'k', 9);
%Apply the Hermite reconstruction for sparse data
result_her11 = hermite_reconstruction(data, dt, 'prev', legpoints, 'solver', 'fmincon', 'search_agents', 5);
mu = @(x,par)par(1).*x;
sigma = @(x,par)par(2)+0 .*x; %this is true model  
par = zeros(2, 1);
par(1) = -1;
par(2) = 1; %true model parameters
xplot = linspace(L, R, 2000); % a dense mesh across the considered range
plot_results(result11, xplot, mu(xplot, par), sigma(xplot, par)); % Euler & true models 
plot_results(result_her11, xplot, mu(xplot, par), sigma(xplot, par)); %Hermite & true models 

