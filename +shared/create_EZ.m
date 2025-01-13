function EZ1 = create_EZ(J, K, ndiff_mu, ndiff_sigma, filename, funname)
    %this function can create a matlab function that solves the Hermite+Taylor expansion
    %of the refined Hermite method (refined Aït-Sahalia method)
    %
    %see equation 10 (appendix):
    %E(Z) = ...
    %J is the order of the Hermite expansion
    %K is the order of the Taylor expansion
    %ndiff_mu are the number of non-zero derivates of mu (can be Inf)
    %ndiff_sigma are the number of non-zero derivatives of sigma (can be Inf)
    %
    %The function returns a function handle. This function takes dt and the mu 
    %and sigma+plus all relevant derivates
    if nargin < 5
        filename = 'EZ';
    end
    if nargin < 6
        funname = filename;
    end
    %comment for the resulting function
    comm=sprintf('E(Z) function of refined Hermite, J=%d  K=%d ndiff_mu = %d, ndiff_sigma = %d',J,K,ndiff_mu,ndiff_sigma);
    
    fprintf('Creating new EZ function, can take some time\n%s\n',comm);
    
    %create symbolic vector
    EZ = sym('EZ', [J 1]);
    syms u m(u) s(u) x x0 dt Ds0 Dm0
 
    zerodiffs = sym([]);

    if ndiff_sigma == 0
        %shortcut if we have additive noise, we replace the function by a
        %constant
        s(u) = Ds0;
    elseif ~isinf(ndiff_sigma)
        %ndiff_sigma is the number of non-negative derivatives, we replace
        %the
        ds = diff(s(x), ndiff_sigma + 1); %#ok<NODEF>
        zerodiffs = [zerodiffs ds];% diff(ds, x)];
    else
        ndiff_sigma = 2 * K - 2;
    end
    if ndiff_mu == 0
        m(u) = Dm0;
    elseif ~isinf(ndiff_mu)
        dm = diff(m(x), ndiff_mu + 1); %#ok<NODEF>
        zerodiffs = [zerodiffs dm];% diff(dm, x)];
    else
        ndiff_mu = 2 * K - 2;
    end
    zero_vals = zeros(size(zerodiffs));
    B = int(1 / s, u, x0, x);
    %for additive noise this equals B=(x - x0)/s(x)
 
    %***************************Main loop ************************
    for j = 1:J
        B0 = B^j;
        S = B0;
        for k = 1:K
            B0 = diff(B0, x);
            %substitute zero differentials with zero to simplify equation
            B0 = subs(B0, zerodiffs, zero_vals);
            B0 = simplifyFraction(m(x) * B0 + 1 / 2 * s(x)^2 * diff(B0, x));
            S = S + B0 * dt^k / shared.factorial1(k);
        end
        EZ(j) = dt^(-j / 2) * limit(S, x, x0);
    end
    EZ = simplifyFraction(EZ);
    
    %replace differentials with variables
    if ndiff_mu > 0
        %usually ndiff_mu>0
        Dm = sym('Dm', [1 ndiff_mu]);
        for i = ndiff_mu: -1:1
            EZ = subs(EZ, {diff(m(x0), i)}, {Dm(i)});
        end
        syms Dm0
        EZ = subs(EZ, m(x0), Dm0);
    else
        Dm = sym([]);
    end
    if ndiff_sigma > 0
        %multiplicative noise
        Ds = sym('Ds', [1 ndiff_sigma]);

        for i = ndiff_sigma: -1:1
            EZ = subs(EZ, {diff(s(x0), i)}, {Ds(i)});
        end
        %in case of additive noise Ds0 is already in EZ
        syms Ds0
        EZ = subs(EZ, s(x0), Ds0);
    else
        Ds = sym([]);
    end
    
    %some final optimizations and create function
    [n, d] = numden(EZ);
    EZ = horner(n, Ds0) ./ horner(d, Ds0); %based on previous comment you should avoid this line in case of a possible error message
    I = hasSymType(EZ, 'variable');
    EZ(~I) = EZ(~I) + 10^(-100) * Ds0; % this is just to turn constant component of ez into a vector by adding an infinitesimal vector
    Dm = [Dm0 Dm];
    Ds = [Ds0 Ds];
    %The reason we SAVE ez as an m-file in bellow is that this way matlab optimizes ez which is very important to reduce computational burden
    matlabFunction(EZ, 'File', filename, 'Vars', [dt Dm Ds], 'Outputs', {'EZ'},'Comment',comm);
    EZ1= str2func(funname);
end
