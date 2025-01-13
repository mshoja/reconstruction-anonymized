function fun = get_hermite(J, K, nder_mu, nder_sigma, refined)
    if nargin < 5
        refined = true;
    end
    base=mfilename('fullpath');
    f=strfind(base,'\+shared\get_hermite');
    base=base(1:f(1)-1);
    if ~exist([base '\+coefficients'],'dir')
        mkdir([base '\+coefficients']);
        mkdir([base '\+coefficients\+ETA']);
        mkdir([base '\+coefficients\+EZ']);     
    end
    if ~refined
        filename = sprintf('%s\\+coefficients\\+ETA\\ETA_m%d_s%d_%d_%d',base, nder_mu, nder_sigma, J, K);
        funname = sprintf('coefficients.ETA.ETA_m%d_s%d_%d_%d', nder_mu, nder_sigma, J, K);
    else        
        filename = sprintf('%s\\+coefficients\\+EZ\\EZ_m%d_s%d_%d_%d',base, nder_mu, nder_sigma, J, K);
        funname = sprintf('coefficients.EZ.EZ_m%d_s%d_%d_%d', nder_mu, nder_sigma, J, K);
    end
    if exist(filename, 'file')
        fun = str2func(funname);
    else
        if refined
            fun = shared.create_EZ(J, K, nder_mu, nder_sigma, filename, funname);
        else
            fun = shared.create_ETA(J, K, nder_mu, nder_sigma, filename, funname);
    %        error('original method of Aït‐Sahalia not yet implemented');
        end
    end
end
