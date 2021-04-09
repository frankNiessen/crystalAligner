function [state, options,optchanged] = myoutputfcn(options,state,flag)
%GAOUTPUTFCNTEMPLATE Template to write custom OutputFcn for GA.
% [STATE, OPTIONS, OPTCHANGED] = GAOUTPUTFCNTEMPLATE(OPTIONS,STATE,FLAG)
% where OPTIONS is an options structure used by GA.
%
% STATE: A structure containing the following information about the state
% of the optimization:
% Population: Population in the current generation
% Score: Scores of the current population
% Generation: Current generation number
% StartTime: Time when GA started
% StopFlag: String containing the reason for stopping
% Selection: Indices of individuals selected for elite,
% crossover and mutation
% Expectation: Expectation for selection of individuals
% Best: Vector containing the best score in each generation
% LastImprovement: Generation at which the last improvement in
% fitness value occurred
% LastImprovementTime: Time at which last improvement occurred
%
% FLAG: Current state in which OutputFcn is called. Possible values are:
% init: initialization state
% iter: iteration state
% interrupt: intermediate state
% done: final state
%
% STATE: Structure containing information about the state of the
% optimization.
%
% OPTCHANGED: Boolean indicating if the options have changed.
%
%	See also PATTERNSEARCH, GA, GAOPTIMSET
% Copyright 2004-2006 The MathWorks, Inc.
% $Revision: 1.1.6.5 $ $Date: 2007/08/03 21:23:22 $
optchanged = false;
switch flag
    case 'init'
        disp('Starting the algorithm');
        if ~isdir([pwd,'\iterOut']); mkdir([pwd,'\iterOut']); end
    case {'iter','interrupt'}
        disp('Iterating ...')
        fname=[pwd,'\iterOut\',num2str(state.Generation),'.mat'];
        save(fname,'state')
    case 'done'
        disp('Performing final task');
        fname=[pwd,'\iterOut\',num2str(state.Generation),'.mat'];
        save(fname,'state')
end
end
