function [npars, index_mu] = getfundimension(mu, sigma)
    haserror = true;
    nparsmu = 1;
    x = [1; 1];
    while nparsmu < 1000 && haserror
        try
            xx = mu(x, rand(1, nparsmu)); %#ok<NASGU>
            haserror = false;
        catch err
            switch err.identifier
                case 'MATLAB:badsubscript'
                    nparsmu = nparsmu + 1;
                    if nparsmu == 1000
                        error('The function mu is not correct, the arguments should be like this: mu(x,par)\n%s', err.message);
                        %rethrow(err);
                    end
                otherwise
                    rethrow(err)
            end
        end
    end
    npars = nparsmu;
    if isempty(sigma) || (ischar(sigma) && strcmpi(sigma, 'additive'))
        npars = npars + 1;
    else
        haserror = true;
        while npars < 1000 && haserror
            try
                xx = sigma(x, rand(1, npars)); %#ok<NASGU>
                haserror = false;
            catch err
                switch err.identifier
                    case 'MATLAB:badsubscript'
                        npars = npars + 1;
                        if npars == 1000
                            error('The function sigma is not correct, the arguments should be like this: sigma(x,par)\n%s', err.message);
                            %rethrow(err);
                        end
                    otherwise
                        rethrow(err)
                end
            end
        end
    end
    index_mu = false(1, npars);
    index_mu(1:nparsmu) = true;
 
end
