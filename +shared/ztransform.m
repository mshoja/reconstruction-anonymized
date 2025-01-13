function [zdata, knots,options] = ztransform(data, datastats,knots, options)

    if options.zscore
        %normal variables
        zdata = (data - datastats.mean) ./ datastats.std;
        options.l = (options.l - datastats.mean) / datastats.std;
        options.r = (options.r - datastats.mean) / datastats.std;
        if ~isempty(options.knots)
            options.knots = (options.knots - datastats.mean) / datastats.std;
        end
        knots = (knots - datastats.mean) / datastats.std;
        
        %rate variables
        options.lb = options.lb / datastats.std;
        options.ub = options.ub / datastats.std;
    end
end


