function optim = setOptimOpts(optim)
fprintf('-> Setting up optimization algorithm ...\n');                 
if optim.order >= 2
    optim.Alg = 'gamultiobj';                                              %'gamultiobj' - multiobjective optimization
elseif optim.order == 1
    optim.Alg = 'ga';                                                      %'ga' - single objective optimitzation
else
    error('%.0f is not a valid objective number',optim.order);            %Error check
end

%% Checking version of Matlab
v = ver('Matlab');
y = str2double(regexp(v.Release,'\d*','Match'));
%% Using 'optimoptions' or 'gaoptimset' for newer and older Matlab releases
if y >= 2016 %Newer or equal Matlab2016a
    optim.opt = optimoptions(optim.Alg);                                   %Create optimization options
    if optim.order >= 2
        if optim.plot
            optim.opt.PlotFcn = {@gaplotpareto};
        end            %Plot Pareto Front
        if optim.hybridFcn
                optim.opt.HybridFcn = {@fgoalattain};
        end      %Set hybrid function
    elseif optim.order == 1
        if optim.plot
            optim.opt.PlotFcn = {@gaplotbestf};
        end             %Plot best fitness
    end
    if optim.iterOut
         optim.opt.OutputFcn = @myoutputfcn;                                %Outputfunction for saving interation states
    end
else %Older Matlab version: use gaoptimset
    optim.opt = gaoptimset(optim.Alg);                                     %Create optimization options
     if optim.order >= 2
        if optim.plot
            optim.opt.PlotFcns = {@gaplotpareto};
        end            %Plot Pareto Front
        if optim.hybridFcn
                optim.opt.HybridFcns = {@fgoalattain};
        end      %Set hybrid function
    elseif optim.order == 1
        if optim.plot
            optim.opt.PlotFcns = {@gaplotbestf};
        end             %Plot best fitness
     end
    if optim.iterOut
        optim.opt.OutputFcn = @myoutputfcn;                                %Outputfunction for saving interation states
    end
end
%General settings
optim.opt.Display = 'final';                                           %Set to 'none', 'off', 'iter', 'diagnose', or 'final'
optim.opt.PopulationSize = optim.popSz;                                %Set Population size
optim.opt.FunctionTolerance = optim.funcTol;                           %Set Function tolerance
optim.opt.MaxStallGenerations = optim.maxStallGen;                     %Set maximum stall generations
if optim.multiCore
    optim.multiCore = true;
else
    optim.multiCore = false;
end
optim.opt.UseParallel = optim.multiCore;                               %Set parallel processing flag
end