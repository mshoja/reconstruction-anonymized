function results = back_ztransform(results)
    if results.options.zscore
        m = results.datastats.mean;
        s = results.datastats.std;
        %normal variables
        results.options.l = results.options.l * s + m;
        results.options.r = results.options.r * s + m;
        if ~isempty(results.options.knots)
            results.options.knots = results.options.knots * s + m;
        end
        results.knots{1} = results.knots{1} * s + m;
        results.knots{2} = results.knots{2} * s + m;
    
        %rate variables
        results.estimated_par = results.estimated_par * s;
        results.options.lb = results.options.lb * s;
        results.options.ub = results.options.ub * s;
        mufun1=@results.mufun;
        results.mufun= @(x,par)s*mufun1((x-m)/s,par/s);
        sigmafun1=@results.sigmafun;
        results.sigmafun= @(x,par)s*sigmafun1((x-m)/s,par/s);
    end
end
