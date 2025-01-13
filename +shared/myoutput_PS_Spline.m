function stop = myoutput_PS_Spline(varargin)

    optimValues = varargin{1};
    state = varargin{2};
    dt = varargin{3};
    M = varargin{4};
    Max_iter = varargin{5};
    if nargin > 5
        knots = varargin{6};
        mesh_fine = varargin{7};
        ModelType = varargin{8};
        SplineType = varargin{9};
    end

    stop = false;
    par = optimValues.bestx;

    %Converting the estimated parameters to the actual unit of data
    par(1:M) = par(1:M) ./ dt;
    par(M + 1:end) = par(M + 1:end) ./ sqrt(dt);

    if isequal(optimValues.bestfval, realmax)
        fprintf('\nThe algorithm could not find a legitimate solution so far. If this lasts long then you should stop the code, pay attention to the following possible reasons for the failure, and re-run again');
        fprintf('\n1) There is a possibility that you have considered a very big (more probable) or very small (less probable) bound around the Euler solution. You should shrink (or expand) the bounds around the Euler solution');
        fprintf('\n2) There is a possibility that there does not exist enough data near data borders. If so, consider shrinking your domain slightly (do not forget to replace the data which fall outside the new domain by NAN. Be careful not to cut them)');
        fprintf('\n3) If your data does not have a high-resolution then you do not need to use a very small ''K'' value which might');
        fprintf('\n4) Consider using the flag ''GWO'' as your solver which is a stronger solver');
        fprintf('\n   be a possible reason for the failure (even if it works, it is useless for you). So, simply try bigger ''K'' values\n');
        fprintf('\nIf (1),(2),(3),(4) are not the reason for the failure then consider the following actions:');
        fprintf('\n1) Reduce the number of knots (but not so much)');
        fprintf('\n2) Consider a simpler model (e.g., use quadratic splines if you are using cubic splines)');
        fprintf('\n3) A nice option is to use the ''Approximate'' flag as spline type\n\n');
    else
        disp('Estimated parameters : ');
        disp(num2str(par));
        disp(['Approximate value of objective function (negative sum of log-likelihoods): ' num2str(optimValues.bestfval)]);

        if nargin > 5
            close all;
            switch ModelType
                case {'Additive noise1', 'Additive noise2', 'Additive noise'}
                    switch SplineType
                        case 'L'
                            mu = @(x)interp1(knots,par(1:M),x);
                            plot(mesh_fine, mu(mesh_fine), '-k');
                            hold on;
                            plot(mesh_fine, 0 .* mesh_fine, '-k');
                            hold on;
                            plot(knots, par(1:M), '*r');
                        case 'Q'
                            sp_mu = spapi(3, knots, par(1:M));
                            plot(mesh_fine, fnval(sp_mu, mesh_fine), '-k');
                            hold on;
                            plot(mesh_fine, 0 .* mesh_fine, '-k');
                            hold on;
                            plot(knots, par(1:M), '*r');
                        case {'C', 'P', 'SCS', 'Approximate'}
                            switch SplineType
                                case 'C'
                                    pp_mu = spline(knots, par(1:M));
                                case 'P'
                                    pp_mu = pchip(knots, par(1:M));
                                case {'SCS', 'Approximate'}
                                    pp_mu = csaps(knots, par(1:M));
                            end
                            plot(mesh_fine, ppval(pp_mu, mesh_fine), '-k');
                            hold on;
                            plot(mesh_fine, 0 .* mesh_fine, '-k');
                            hold on;
                            plot(knots, par(1:M), '*r');
                            hold on;
                    end
                    annotation('textbox', [0.55, 0.8, 0.1, 0.1], 'String', strcat(['Noise intensity is : ' num2str(par(end))]));
                    xlabel('State, x');
                    ylabel('Drift function, \mu (x)');
                    xlim([-inf inf]);
                    ylim([-inf inf]);
                case 'Multiplicative noise'
                    switch SplineType
                        case 'LL'
                            pp_mu = interp1(knots, par(1:M), 'linear', 'pp');
                            pp_sigma = interp1(knots, par(M + 1:end), 'linear', 'pp');
                            figure(1);
                            plot(mesh_fine, ppval(pp_mu, mesh_fine), '-k');
                            hold on;
                            plot(knots, par(1:M), '*r');
                            hold on;
                            plot(knots, 0 .* knots, '-k');
                            xlabel('State, x');
                            ylabel('Drift function, \mu (x)');
                            xlim([-inf inf]);
                            ylim([-inf inf]);
                            pos1 = get(gcf, 'Position');
                            set(gcf, 'Position', pos1 - [pos1(3) / 2, 0, 0, 0]);
                            figure(2);
                            plot(mesh_fine, ppval(pp_sigma, mesh_fine), '-k');
                            hold on;
                            plot(knots, par(M + 1:end), '*r');
                            xlabel('State, x');
                            ylabel('Diffusion function, \sigma (x)');
                            xlim([-inf inf]);
                            ylim([0 inf]);
                            pos2 = get(gcf, 'Position');
                            set(gcf, 'Position', pos2 + [pos1(3) / 2, 0, 0, 0]);
                            figure(3);
                            A1 = ppval(pp_mu, mesh_fine);
                            A2 = ppval(pp_sigma, mesh_fine);
                            U = -2 .* cumtrapz(mesh_fine, A1 ./ A2.^2) + 2 .* log(A2);
                            U = U - min(U) + 0.5;
                            plot(mesh_fine, U, '-k');
                            xlabel('State, x');
                            ylabel('Effective potential, U_{eff}(x)');
                            xlim([-inf inf]);
                            ylim([0 inf]);
                            pos3 = get(gcf, 'Position');
                            set(gcf, 'Position', pos3 + [0, -504, 4, 0]);
                        case {'CP', 'PP', 'CC', 'SCS', 'Approximate'}
                            switch SplineType
                                case 'CP'
                                    pp_mu = spline(knots, par(1:M));
                                    pp_sigma = pchip(knots, par(M + 1:end));
                                case 'PP'
                                    pp_mu = pchip(knots, par(1:M));
                                    pp_sigma = pchip(knots, par(M + 1:end));
                                case 'CC'
                                    pp_mu = spline(knots, par(1:M));
                                    pp_sigma = spline(knots, par(M + 1:end));
                                case {'SCS', 'Approximate'}
                                    pp_mu = csaps(knots, par(1:M));
                                    pp_sigma = csaps(knots, par(M + 1:end));
                            end
                            figure(1);
                            plot(mesh_fine, ppval(pp_mu, mesh_fine), '-k');
                            hold on;
                            plot(mesh_fine, 0 .* mesh_fine, '-k');
                            hold on;
                            plot(knots, par(1:M), '*r');
                            xlabel('State, x');
                            ylabel('Drift function, \mu (x)');
                            xlim([-inf inf]);
                            ylim([-inf inf]);
                            pos1 = get(gcf, 'Position');
                            set(gcf, 'Position', pos1 - [pos1(3) / 2, 0, 0, 0]);
                            figure(2);
                            plot(mesh_fine, ppval(pp_sigma, mesh_fine), '-k');
                            hold on;
                            plot(knots, par(M + 1:end), '*r');
                            xlabel('State, x');
                            ylabel('Diffusion function, \sigma (x)');
                            xlim([-inf inf]);
                            ylim([0 inf]);
                            pos2 = get(gcf, 'Position');
                            set(gcf, 'Position', pos2 + [pos1(3) / 2, 0, 0, 0]);
                            figure(3);
                            A1 = ppval(pp_mu, mesh_fine);
                            A2 = ppval(pp_sigma, mesh_fine);
                            U = -2 .* cumtrapz(mesh_fine, A1 ./ A2.^2) + 2 .* log(A2);
                            U = U - min(U) + 0.5;
                            plot(mesh_fine, U, '-k');
                            xlabel('State, x');
                            ylabel('Effective potential, U_{eff}(x)');
                            xlim([-inf inf]);
                            ylim([0 inf]);
                            pos3 = get(gcf, 'Position');
                            set(gcf, 'Position', pos3 + [0, -504, 4, 0]);
                        case 'QQ'
                            sp_mu = spapi(3, knots, par(1:M));
                            sp_sigma = spapi(3, knots, par(M + 1:end));
                            figure(1);
                            plot(mesh_fine, fnval(sp_mu, mesh_fine), '-k');
                            hold on;
                            plot(mesh_fine, 0 .* mesh_fine, '-k');
                            hold on;
                            plot(knots, par(1:M), '*r');
                            xlabel('State, x');
                            ylabel('Drift function, \mu (x)');
                            xlim([-inf inf]);
                            ylim([-inf inf]);
                            pos1 = get(gcf, 'Position');
                            set(gcf, 'Position', pos1 - [pos1(3) / 2, 0, 0, 0]);
                            figure(2);
                            plot(mesh_fine, fnval(sp_sigma, mesh_fine), '-k');
                            hold on;
                            plot(knots, par(M + 1:end), '*r');
                            xlabel('State, x');
                            ylabel('Diffusion function, \sigma (x)');
                            xlim([-inf inf]);
                            ylim([0 inf]);
                            pos2 = get(gcf, 'Position');
                            set(gcf, 'Position', pos2 + [pos1(3) / 2, 0, 0, 0]);
                            figure(3);
                            A1 = fnval(sp_mu, mesh_fine);
                            A2 = fnval(sp_sigma, mesh_fine);
                            U = -2 .* cumtrapz(mesh_fine, A1 ./ A2.^2) + 2 .* log(A2);
                            U = U - min(U) + 0.5;
                            plot(mesh_fine, U, '-k');
                            xlabel('State, x');
                            ylabel('Effective potential, U_{eff}(x)');
                            xlim([-inf inf]);
                            ylim([0 inf]);
                            pos3 = get(gcf, 'Position');
                            set(gcf, 'Position', pos3 + [0, -505, 4, 0]);
                    end
            end
            drawnow;
        end
    end
end















