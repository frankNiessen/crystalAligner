function [oNew,stgRot,x,eps] = crystalAligner(crys,stg,optim,FIB)
%[oNew,stgRot,x,eps] = crystalAligner(crys,stg,optim,FIB)
% crystalAligner main function
% This function is generally called from an example script
% See example scripts and README.md for information on input parameters

%% Setup Optimization options
crys = checkerror(crys);                                                   %Error checking
optim = setOptimOpts(optim);                                               %Optimization initialization function
scrPrnt('Ini',crys,stg,optim);                                             %Screen print of optimization objectives and limits parameters

%% Plot initial orientations
plotOrientations(optim,{crys.ori},crys,'initial');

%% Optimization - Multiobjective genetic algorithm
[oNew,stgRot,x,eps] = runOptim(crys,stg,optim,FIB);                        %Optimization function

%% Plot stereographic projection and inverse polefigure of aligned equivalent crystal directions
plotOrientations(optim,oNew,crys,'result');                                %Plotting the results
end

