function results = simlangevin(mu, sigma, pars, ts, x0, deltat, limits)
    %This is the main trick of the fitlangevin method. Here an efficient vectorized run
    %is done, x0 is typically a long vector that is run simultaneously.
    % 
    if nargin < 4 || isempty(ts)
        ts = (1:1000)';
    end
    if nargin < 5 || isempty(x0)
        x0 = 0.001;
    end
    if nargin < 6 || isempty(deltat)
        deltat = 0.01;
    end
    if nargin < 7
        limits = [];
    end
    if all(isinf(limits))
        limits = [];
    end

    x = x0(:)';
    sqrt_deltat = sqrt(deltat);
    t = ts(1):deltat:ts(end) + deltat;
    results = zeros(numel(ts), numel(x));
    results(1, :) = x;
    siz = size(x);
    nextt = 2;
    for i = 2:numel(t)
        x = x + mu(x, pars) .* deltat + sigma(x, pars) .* randn(siz) .* sqrt_deltat;
        if ~isempty(limits)
            x(x < limits(1)) = limits(1);
            x(x > limits(2)) = limits(2);
        end
        if nextt <= numel(ts) && t(i) >= ts(nextt)
            results(nextt, :) = x;
            nextt = nextt + 1;
        end
    end
end

