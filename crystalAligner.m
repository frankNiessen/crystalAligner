%% crystalAligner.m - Crystal Alignment Tool for Electron Microscopy 
% *************************************************************************
% Script for optimized alignment of one or two crystal directions of
% given crystal orientations with major coordinate axes of a scanning 
% electron microscope 
% Please find more information and examples in the attached documentation
% *************************************************************************
% Dr. Frank Niessen, University of Wollongong, Australia, 2018
% contactnospam@fniessen.com (remove the nospam to make this email address 
% work)
% Published on 11/07/2019 - version 1.0
% *************************************************************************
% Execution of the software requires installation of 
% * mtex (http://mtex-toolbox.github.io/)
% * MATLAB Optimization toolbox
% * MATLAB Global optimization toolbox
% * Optional: MATLAB Parallel Computing Toolbox
% Compatibility tested with MATLAB R2104a - R2016b and MTEX 4.2.1 to 5.1.1 
% *************************************************************************
% Copyright 2018 Frank Niessen (see attached MIT license file)

function [epsilon,rot,oNew] = crystalAligner()
%function crystalAligner()
%epsilon:   Angular deviation from ideal alignment (fitness function value)
%rot:       Stage rotation angles of optimal solution (optimal individuals of fitness function
%oNew:      New orientation after alignment of crystal
clc; close all; warning('off','all');
fprintf('\n*************************************************************');
fprintf('\n                      crystalAligner v. 1.0 \n');
fprintf('*************************************************************\n\n');
fprintf('\n -> Starting up MTEX ...');                                     %ScreenPrint
iniExtLibs;                                                                %Automatically open and initialize MTEX and check for MATLAB toolboxes
%**************************************************************************
%% Initialization
% ****************************** Crystals *********************************
% *** Crystal Alignment Objective 1
crys.o(1,:)     = [61 42 9];                                               %Crystal orientation in Euler angles [pih1 Phi phi2]       
crys.cs{1}      = 'cubic';                                                 %Crystal structure string (follow MTEX convention)
crys.alignAx(1) = yvector;                                                 %Microscope axis for alignment with crystal direction/plane; Examples: zvector; [.5 .5 1]; xvector; ...
crys.Miller{1}  = [1 1 3];                                                 %Miller indices for alignment (in Multiobjective Optimization several Miller-sets will start several optimizations);  Examples: [1 0 0; 1 1 0]; [-1 2 1]; ...
crys.type{1}    = 'uvw';                                                   %Type of Miller: 'hkl': Crystal plane; 'uvw': Crystal direction
crys.sym(1)     = 1;                                                       %Apply crystal symmetry: 1: yes 0: no
% *** Crystal Alignment Objective 2
crys.o(2,:)     = [61 42 9];                                               %Crystal orientation in Euler angles [pih1 Phi phi2]       
crys.cs{2}      = 'cubic';                                                 %Crystal structure string (follow MTEX convention)
crys.alignAx(2) = zvector;                                                 %Microscope axis for alignment with crystal direction/plane; Examples: zvector; [.5 .5 1]; xvector; ...
crys.Miller{2}  = [1 1 0];                                                 %Miller indices for alignment (in Multiobjective Optimization several Miller-sets will start several optimizations);  Examples: [1 0 0; 1 1 0]; [-1 2 1]; ...
crys.type{2}    = 'uvw';                                                   %Type of Miller: 'hkl': Crystal plane; 'uvw': Crystal direction
crys.sym(2)     = 1;                                                       %Apply crystal symmetry: 1: yes 0: no
% ******************************* Stage ***********************************                                          
stg.rot     = [xvector; zvector];                                          %Stage rotation axes                                        
stg.LB      = [    0     -180  ];                                          %Lower bound [°]
stg.UB      = [   20      180  ];                                          %Upper bound [°]                                         
stg.sign    = [    1      -1   ];                                          %Sign 1: Right hand rule convention; -1: Left hand rule convention
stg.order   = [    2      1    ];                                          %Hierarchy / order of rotation: Apply rotation around axis 2 and then 1; Example: [2 1];
% ************************* Genetic algorithm *****************************
%Genetic algorithm
optim.popSz = 100;                                                         %Population size
optim.funcTol = 0.1;                                                       %FunctionTolerance
optim.maxStallGen = 10;                                                    %Maximum stall generations
optim.iterOut = 0;                                                         %Writing output for each iteration in subFolder 'iterOut'
%Multiobjective genetic algorithm settings
optim.wghtFac = [1,1];                                                     %Weighting factors for TOPSIS multiobjective decision making method
optim.multiCore = 0;                                                       %Flag: Utilization of parallel processing (switch off if errors ocur) [1|0]
optim.hybridFcn = 0;                                                       %Flag: Use a hybrid function to (may speed op convergence but compromise diversity of solution space) [1|0]
optim.autoSol = 1;                                                         %Flag: Pick optimum solution automatically by distance of Pareto solution from the optimal solution [1|0]
% ***************************** Optional **********************************
%FIB liftout calculations
FIB.mode = 1;                                                              %Flag: FIB liftout output [1|0]
FIB.trench.ang = 52;                                                       %Trenching - or look-in - angle of Trench [°]
FIB.trench.z = 15;                                                         %Trench depth 'z' [µm]
FIB.axs.tilt = 1;                                                          %Index of tilt axis in 'axs.rot'
FIB.axs.rot = 2;                                                           %Index of rotation axis in 'axs.rot'                                                        
%Output
optim.plot = 1;                                                            %Plotting 1: On 0: Off
%
%
%
%
%% Error checking and PreProcessing
crys = checkerror(crys);                                                   %Error checking
%PreProcessing
crys.nrObj = length(crys.Miller);                                          %Number of alignment objectives
crys.ss = 'orthorhombic';                                                  %Specimen system
crys.MillerAll = unique(vertcat(crys.Miller{:}),'rows');                   %All Miller indicees
%% Setup Optimization options
optim = setOptimOpts(crys,optim);                                          %Optimization initialization function
scrPrnt('Ini',crys,stg,optim);                                             %Screen print of optimization objectives and limits parameters
%% Define Crystal system, orientation and crystal directions in mtex & Plot orientations
[o,cs,ss,crys] = defOri(crys,optim);                                       %Orientation definition function
%% Optimization - Multiobjective genetic algorithm 
[oNew,stgRot,x,eps] = runOptim(crys,stg,optim,o,FIB,ss);                    %Optimization function
%% Plot stereographic projection and inverse polefigure of aligned equivalent crystal directions
plotting(optim,oNew,crys);                                                 %Plotting function
epsilon = eps.opt;                                                         %Assign return value
rot = x.opt;                                                               %Assign return value
fprintf('\n -> All done!\n');                                              %ScreenPrint
end
%% iniMtex
function iniExtLibs
%Initialize MTEX
try
    MTEXmenu;                                                              %Open m-tex menu
catch
    startup_mtex;                                                          %Startup m-tex
end   
% Check for toolboxes
assert(logical(license('test','gads_toolbox')),'Global Optimization Toolbox not installed'); 
assert(logical(license('test','optimization_toolbox')),'Optimization Toolbox not installed');
end
%% checkError - Error checking of input data
function crys = checkerror(crys)
assert(size(crys.alignAx,2) == length(crys.Miller), 'Number of axis-crystalDirection pairs of alignment objetives not equal');  %Check number of alignment axes and alignment objectives
assert(length(crys.cs) == size(crys.o,1),'Number of CrystalSystem-Orientation pairs not equal');  %Check number of crystal systems and orientations
%Check for different Miller notations
for i = 1:length(crys.Miller) 
   assert(any(strcmpi(crys.type{i},{'hkl','hkil','uvw','uvtw'})),'Please choose ''uvw'', ''uvtw'', ''hkl'' or ''hkil'' as valid Miller type');
   assert(size(crys.Miller{i},2) == 3 || size(crys.Miller{i},2) == 4,'Only 3 or 4 Miller indicees per set allowed');
   if size(crys.Miller{i},2) == 4 && any(strcmpi(crys.type(i),{'hkl','hkil'})) %hkil notation
      crys.Miller{i} = crys.Miller{i}(:,[1 2 4]);
      crys.type{i} = 'hkl';
   elseif size(crys.Miller{i},2) == 4 && any(strcmpi(crys.type(i),{'uvw','uvtw'})) %uvtw notation
      crys.Miller{i} = [crys.Miller{i}(:,1)-crys.Miller{i}(:,3),crys.Miller{i}(:,2)-crys.Miller{i}(:,3),crys.Miller{i}(:,4)];   
      crys.type{i} = 'uvw';
   end    
end
end
%% setOptimOpts - Setting optimization solver options
function optim = setOptimOpts(crys,optim)
fprintf('\n -> Setting up optimization algorithm ...');                    %ScreenPrint
if crys.nrObj >= 2 
    optim.Alg = 'gamultiobj';                              %'gamultiobj' - multiobjective optimization
elseif crys.nrObj == 1 
    optim.Alg = 'ga';                                  %'ga' - single objective optimitzation
else
    error('%.0f is not a valid objective number',crys.nrObj);             %Error check
end

%% Checking version of Matlab
v = ver('Matlab');
y = str2double(regexp(v.Release,'\d*','Match'));
%% Using 'optimoptions' or 'gaoptimset' for newer and older Matlab releases
if y >= 2016 %Newer or equal Matlab2016a
    optim.opt = optimoptions(optim.Alg);                                   %Create optimization options   
    if crys.nrObj >= 2
        if optim.plot
            optim.opt.PlotFcn = {@gaplotpareto}; 
        end            %Plot Pareto Front
        if optim.hybridFcn
                optim.opt.HybridFcn = {@fgoalattain}; 
        end      %Set hybrid function
    elseif crys.nrObj == 1
        if optim.plot 
            optim.opt.PlotFcn = {@gaplotbestf}; 
        end             %Plot best fitness   
    end
    if optim.iterOut
         optim.opt.OutputFcn = @myoutputfcn;                                %Outputfunction for saving interation states
    end
else %Older Matlab version: use gaoptimset
    optim.opt = gaoptimset(optim.Alg);                                     %Create optimization options
     if crys.nrObj >= 2
        if optim.plot
            optim.opt.PlotFcns = {@gaplotpareto}; 
        end            %Plot Pareto Front
        if optim.hybridFcn
                optim.opt.HybridFcns = {@fgoalattain}; 
        end      %Set hybrid function
    elseif crys.nrObj == 1
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
optim.opt.UseParallel = optim.multiCore;                               %Set parallel processing flag
end

%% function defOri - Defining and plotting of orientations
function [o,cs,ss,crys] = defOri(crys,optim)
fprintf('\n -> Defining crystal system, orientation and directions ...');  %ScreenPrint
ss = crystalSymmetry(crys.ss);                                             %Define Specimen system
for c = 1:size(crys.o,1) %Loop over crystals
    cs{c} = crystalSymmetry(crys.cs{c});                                   %Define Crystal system    
    o{c} =  orientation('Euler',crys.o(c,1)*degree,crys.o(c,2)*degree,...
                                crys.o(c,3)*degree,cs{c},ss);              %Compute crystal orientation
    %Define Miller planes or directions
    for i = 1:size(crys.Miller{c},1) %Loop over Miller indicee sets
        crys.oMil{c}{i} = Miller(crys.Miller{c}(i,1),crys.Miller{c}(i,2),...
                                 crys.Miller{c}(i,3),cs{c},crys.type{c});  %Create Miller indexed direction
        crys.strMil{c}{i} = regexprep(num2str(round(...
                                      crys.oMil{c}{i}.(crys.type{c}))),'\s+',''); %Save Miller labels
    end 
end

% *** Plot stereographic projection and inverse polefigure of equivalent crystal directions
% Reference frame is stage coordinate system C_s or equivalently microscope
% coordinate system C_m for 0° stage tilt, no pre-tilt and no change in
% rotation with respect to the orientation measurement
if optim.plot
    for c = 1:size(crys.o,1)
        fprintf(['\n -> Plotting stereographic projections and IPF of',...
                 ' crystal %.0f ...'],c);                                   %ScreenPrint
        stereoProj(o{c},crys.oMil{c},crys.strMil{c},crys.sym(c),...
                   sprintf('Original orientation %.0f',c));                %Plot stereographic projection of equivalent crystal directions of crystal orientation 'o'
        figure;                                                            %Create figure
        plotIPDF(o{c},[xvector,yvector,zvector],'antipodal',...
                 'MarkerFaceColor','k');                                   %Plot inverse pole-figure (IPF) [x:(100), y:(010) and z:(001)]
        set(gcf,'name',sprintf('Original orientation %.0f',c));
    end
tileFigs();                                                                %Sort figures
end
end

%% runOptim - Optimization function
function [oNew,stgRot,x,eps] = runOptim(crys,stg,optim,o,FIB,ss)
for i = 1:length(crys.oMil{1}) %Loop over crystal directions
    scrPrnt('Optim',crys,optim,i);                                     %Screen output optimization problem
    % *** Multiobjective opimization **************************************
    if crys.sym(1)
        dirs = symmetrise(crys.oMil{1}{i}); 
    else
        dirs = crys.oMil{1}{i};
    end
    if strcmp(optim.Alg,'gamultiobj')
        % *** Define optimization function *****
        fMin = @(x)minFunc2D(x,o,dirs,crys.oMil{2},stg,crys);
        % *** Run optimization *****
        % Checking version of Matlab
        v = ver('Matlab');
        y = str2double(regexp(v.Release,'\d*','Match'));
        %% Using 'optimoptions' or 'gaoptimset' for newer and older Matlab releases
        if y >= 2015 %Newer or equal Matlab2016a
            [x.ga,eps.ga] = gamultiobj(fMin,size(stg.rot,1),[],[],[],[],stg.LB,stg.UB,[],optim.opt);
        else
            [x.ga,eps.ga] = gamultiobj(fMin,size(stg.rot,1),[],[],[],[],stg.LB,stg.UB,optim.opt);
        end
        % Choose optimal Pareto solution
        if optim.autoSol
            scrPrnt('Solution',optim,crys);                                 %Screen Output
            ind = topsis(eps.ga,optim.wghtFac,'matrix');                   %TOPSIS multi-objective decision-making
            ind = ind(1);                                                  %Take best solution
        else
             ind = listdlg('PromptString',...
                    'Choose a Pareto-solution: ax{1}[°] ax{2}[°] eps1[°] eps2[°]:',...
                     'liststring',num2str([x.ga,eps.ga]),...
                     'SelectionMode','single','ListSize',[300,250]);       %Choose Pareto-solution manually
        end
        x.opt(i,:) = x.ga(ind,:);                                          %Save optimal solution
        eps.opt(i,:) = eps.ga(ind,:);                                      %Save misalignment of optimal solution
        x.out(i,:) = x.opt(i,:).*stg.sign;                                 %Adapt possible different stage rotation convention for output  

        %Plot optimal solution
        if optim.plot
            h.ax = findobj('type','axes','tag','gaplotpareto');            %Find axes
            hold(h.ax,'on');                                               %Hold on                                          
            plot([0,eps.opt(i,1)],[0,eps.opt(i,2)],'k');                   %Plot line
            plot(eps.opt(i,1),eps.opt(i,2),'*k');                          %Plot marker
        end
        %Find out which second crystal direction was aligned      
        rotTot = rotation('Euler',0,0,0);
        for r = 1:size(stg.order,2)
            axNr = stg.order(r);
            rot(axNr) = rotation('axis',stg.rot(axNr),'angle',x.opt(axNr)*degree);     %Compute rotation around microscope rotation axis 'r' 'stg.rot(r)'
            rotTot = rot(axNr)*rotTot;
        end       
        for j=1:length(crys.oMil{2})
            if crys.sym(2)
                dirs = symmetrise(crys.oMil{2}{j});
            else
                dirs = crys.oMil{2}{j};
            end  
           epsTemp(j)=min(angle(rotTot*o{2}*dirs,crys.alignAx(2))/degree); %Find misalignment of microscope axis 2 'crys.alignAx(2)' with 'o*CD' subject to rotations 'rotX*rotZ'
        end
        [~,crystDirInd] = min(epsTemp);
        crys.str.optAx = crys.strMil{2}{crystDirInd}; 
    % *** Singleobjective opimization *************************************
    elseif strcmp(optim.Alg,'ga')
        % *** Define optimization function *****
        fMin = @(x)minFunc(x,o,dirs,stg,crys);           
        % *** Run optimization *****
        [x.opt(i,:),eps.opt(i,:)] = ga(fMin,size(stg.rot,1),[],[],[],[],stg.LB,stg.UB,[],optim.opt); %genetic algorithm
         x.out(i,:) = x.opt(i,:).*stg.sign;                                %Adapt possible different stage rotation convention for output  

    else
        error('Invalid choice of optimization algorithm');                 %No valid choice of alcrys.oithm
    end
                           %Optimized crystal direction with second axis
    %Output
    scrPrnt('Result',dir,stg,x,eps,i,crys,optim);                          %General results screen output
    scrPrnt('FIB_FEIHelios',x,i,FIB);                                      %Instrument specific screen output
    rotTot = rotation('Euler',0,0,0);                                                            %Ini
    for r = 1:size(stg.rot,1)
        axNr = stg.order(r);
        stgRot{i}.ax{axNr} = rotation('axis',stg.rot(axNr),'angle',x.opt(i,axNr)*degree);  %Convert rotation around microscoe axis r 'stg.rot(r)' to Euler angles
        rotTot = stgRot{i}.ax{axNr}*rotTot;                                   %Total rotation
    end
    for c = 1:length(o)
        oNew{i}{c} = rotTot*o{c};                                          %Compute new crystal orientation after applied stage tilt and rotation
        oNew{i}{c}.SS = ss;                                                %Assign specimen symmetry
    end
end
end
%% minFunc2D - 2D Minimization objective function
function eps = minFunc2D(x,o,CD1,CD2,stg,crys)
%function err = minFunc2D(x,o,CD,tensD,wghtFac)
%Objective function calculating the misalignment of the microscope  
%axes 1 and 2 in 'crys.alignAx' with crystal directions 'CD1' and 'CD2'
%rotated by tilt/rotation angles around microscope axes 1 and 2 in 'stg.rot' 
%under consideration of weighting factors 'wFac'
warning('off','all');
rotTot = rotation('Euler',0,0,0);
for r = 1:size(stg.order,2)
    axNr = stg.order(r);
    rot(axNr) = rotation('axis',stg.rot(axNr),'angle',x(axNr)*degree);     %Compute rotation around microscope rotation axis 'r' 'stg.rot(r)'
    rotTot = rot(axNr)*rotTot;
end
eps(1) = min(angle(rotTot*o{1}*CD1,crys.alignAx(1))/degree);                  %Find misalignment of microscope axis 1 'crys.alignAx(1)' with 'o*CD' subject to rotation 'rotTot'
epsTemp = zeros(length(CD2),1);                                            %Preallocate memory for 'epsTemp'
for i=1:length(CD2)
    if crys.sym(2)
        dirs = symmetrise(CD2{i});
    else
        dirs = CD2{i};
    end    
   epsTemp(i)=min(angle(rotTot*o{2}*dirs,crys.alignAx(2))/degree);     %Find misalignment of microscope axis 2 'crys.alignAx(2)' with 'o*CD' subject to rotations 'rotX*rotZ'
end
eps(2)=min(epsTemp);                                                       %Take the smallest misalignment in axis 2 as optimum solution
end
%% minFunc - 1D Minimization objective function
function err = minFunc(x,o,CD,stg,crys)
%function err = minFunc(x,o,CD)
%Objective function calculating the misalignment of microscope 
%axis 1 'crys.alignAx(1)' with crystal directions 'CD' rotated by 
%tilt/rotation angles around microscope axes 1 and 2 in 'stg.rot'
%x: Rotation angles
%o: Crystal orientation
%CD: Crystal directions
rotTot = rotation('Euler',0,0,0);
for r = 1:size(stg.order,2)
    axNr = stg.order(r);
    rot(axNr) = rotation('axis',stg.rot(axNr),'angle',x(axNr)*degree);     %Compute rotation around microscope rotation axis 'r' 'stg.rot(r)'
    rotTot = rot(axNr)*rotTot;
end
err = min(angle(rotTot*o{1}*CD,crys.alignAx(1))/degree);                      %Find misalignment of microscope axis 1 'crys.alignAx(1)' with 'o*CD' subject to rotation 'rotTot'
end
%% TOPSIS - Technique for Order of Preference by Similarity to Ideal Solution (TOPSIS)
function ind = topsis(epsIn,w,normMode)
%function ind = topsis(epsIn,w)
%TOPSIS - Technique for Order of Preference by Similarity to Ideal Solution (TOPSIS)
%Multi-objective decision-making method to rank Pareto Solutions
%INPUT
%epsIn: Fitness function values [NrParetoSolutions x NrObjectives]
%w:     Objective weights [1 x NrObjectives]
%ind:   Indices of Pareto solutions sorted from highest to lowest score
if strcmpi(normMode,'Vector')
   epsNorm = epsIn./sqrt(sum(epsIn.^2));                                   %Normalized fitness function values by column
elseif strcmpi(normMode,'Matrix') 
   epsNorm = epsIn./sqrt(sum(epsIn(:).^2));                                %Normalized fitness function values by matrix
else
   error(sprintf('Invalid normalization mode ''%s''',normMode));           %Error
end
                                       
M = epsNorm.*repmat(w,size(epsNorm,1),1);                                  %Compute decision matrix                                               
Mworst = max(M,[],1);                                                      %Find maxima
Mbest = min(M,[],1);                                                       %Find minima
Lworst = sqrt(sum((M-repmat(Mworst,size(M,1),1))'.^2));                    %Compute L-square-distance between M and Mworst
Lbest = sqrt(sum((M-repmat(Mbest,size(M,1),1))'.^2));                      %Compute L-square-distance between M and Mbest
score = Lworst./(Lworst+Lbest);                                            %Compute score
[~,ind] = sort(score,'descend');                                           %Rank score
end

%% plotting - Plotting function
function plotting(optim,oNew,crys)
if optim.plot   
    for c = 1:size(crys.o,1)
        for i = 1:length(crys.oMil{c})
            titleStr = sprintf('%s aligned with %s %s',xyzStr(crys.alignAx(c)), crys.oMil{c}{1}.CS.LaueName, crys.strMil{c}{i});
            stereoProj(oNew{i}{c},crys.oMil{c},crys.strMil{c},crys.sym(c),titleStr);%Plot new crystal directions of interest for z-axis
            figure;                                                        %Create figure
            plotIPDF(oNew{i}{c},[xvector,yvector,zvector],'antipodal',...
                                                       'MarkerFaceColor','k'); %Plot inverse pole-figure (IPF) [x:(100), y:(010) and z:(001)]
            set(gcf,'name',titleStr);                                      %set title
        end
    end
    tileFigs();                                                                %Align figures
end
end
%% stereoProj - Plotting of stereographic projection
function stereoProj(o,Dir,Dirstrs,sym,titleStr)
%function stereoProj(o,Dir,Dirstrs,titleStr)
%Plot stereographic projection of all equivalent crystal directions to 
%Miller crystal direction 'Dir' with respect crystal orientation 'o'.
%'Dirstrs' are legend strings and 'titleStr' defines the title string
%% Initialization
figure;                                                                    %Create new figure
markerSz = 12;                                                             %Define marker size
markers = {'o','d','s','v','h','^'};                                       %Define marker style
colors = {'k','r','b','g','m','c'};                                        %Define marker colors
lgFntSz = 20;                                                              %FontSize for legend
%% Plotting
for i = 1:length(Dir) %Loop over crystal directions
    if sym
        dirs = symmetrise(Dir{i});
    else
        dirs = Dir{i};
    end  
    plot(o*dirs,'antipodal','grid','backgroundcolor','w',...
                 'MarkerSize',markerSz,'MarkerEdgeColor','k',...
                 'MarkerFaceColor',colors{i},'Marker',markers{i});         %Plot equivalent crystal directions    
    hold on                                                                %Allow additive plotting
    annotate([xvector,yvector],'label',{'X','Y'},'backgroundcolor','w');   %Annotate labels
    L(i) = plot(nan(6,1),nan(6,1),'MarkerSize',markerSz,...
           'MarkerEdgeColor','k','MarkerFaceColor',colors{i},...
           'Marker',markers{i},'Color','w');                               %Create 'Ghost plot' for creating a nice legend
end
%% Finishing
l = legend(L,Dirstrs,'FontName','TimesNewRoman','FontSize',20);            %Create legend
l.FontSize = lgFntSz;                                                      %Change legend fontsize
set(gcf,'name',titleStr);
end
%% scrPrnt - Screen print of optimization problem
function scrPrnt(mode,varargin)
%function scrPrnt(crys,stg,dir,optim,mode)
switch mode
    case 'Ini'
        %Variable Input
        crys = varargin{1}; 
        stg = varargin{2}; 
        optim = varargin{3};
        fprintf('\n-------------------------------------------------------------\n');
        fprintf('*** Input parameters ***\n');
        fprintf('\nOptimization with algorithm ''%s''',optim.Alg);
        fprintf('\nParallel alignment of:\n');
        for i = 1:size(crys.alignAx,2)
           if crys.sym(i); sym=''; else; sym = '(NoSymmetry)'; end
           fprintf(' - Microscope %s with %s crystal %s: \t%s %s\n',...
                xyzStr(crys.alignAx(i)),crys.cs{i},crys.type{i},num2str(crys.Miller{i}'),sym);            
        end
        fprintf('Rotational degrees of freedom:\n');
        for r = 1:size(stg.rot,1)
            fprintf(' - Rotation around microscope %s: %.1f to %.1f °\n',...
                    xyzStr(stg.rot(r)),stg.LB(r),stg.UB(r));
        end
        fprintf('-------------------------------------------------------------');
    case 'Optim'
        %Variable Input
        crys = varargin{1};
        optim = varargin{2};
        i = varargin{3};
        %Screen output
        fprintf('\n*************************************************************');
        fprintf('\n Optimization problem %.0f - Parallel alignment of:\n',i);
        fprintf('   -> %s %s %s with microscope %s\n',crys.cs{1},crys.strMil{1}{i},crys.type{1},xyzStr(crys.alignAx(1)));
        if strcmp(optim.Alg,'gamultiobj')
           if size(crys.Miller{2},1)==1
               fprintf('   -> %s %s %s with microscope %s\n',crys.cs{2},num2str(crys.Miller{2}'),crys.type{2},xyzStr(crys.alignAx(2))); 
           else
               fprintf('   -> Either of %s %ss with microscope %s\n',crys.oMil{2}{1}.CS.LaueName,crys.type{2}, num2str(crys.Miller{2}'),xyzStr(crys.alignAx(2)));  
           end
        end
        fprintf('*************************************************************\n\n');
    case 'Solution'
        %Variable Input
        optim =  varargin{1}; 
        crys =  varargin{2}; 
        %Screen output
        fprintf('\n*************************************************************');
        fprintf('\nOptimal solution found by shortest distance to origin \n');
        if strcmp(optim.Alg,'gamultiobj')
            fprintf('under consideration of weighting factors %.0f and %.0f for alignment with microscope %s and %s\n',optim.wghtFac(1),optim.wghtFac(2),xyzStr(crys.alignAx(1)),xyzStr(crys.alignAx(2)));
        end
   case 'Result'
        %Variable input
        dir = varargin{1}; 
        stg = varargin{2}; 
        x = varargin{3}; 
        eps = varargin{4};
        i = varargin{5};
        crys = varargin{6};
        optim = varargin{7};      
        %Screen Output
        fprintf('\n-------------------------------------------------------------\n');
        fprintf('*** Optimization results ***');
        if strcmp(optim.Alg,'ga')
            fprintf('\nOptimal parallel alignment of\n\t-microscope %s with crystal %s %s:\n',...
                     xyzStr(crys.alignAx(1)),crys.type{1},crys.strMil{1}{i}); 
        elseif strcmp(optim.Alg,'gamultiobj')
            fprintf('\nOptimal parallel alignment of\n\t- microscope %s with %s crystal %s %s and \n\t- microscope %s with %s crystal %s %s:\n',...
                     xyzStr(crys.alignAx(1)),crys.cs{1},crys.type{1},crys.strMil{1}{i},xyzStr(crys.alignAx(2)),crys.cs{2},crys.type{2},crys.str.optAx); 
        end
        
        for r = 1:size(stg.rot,1)
            fprintf('   -> Rotation around microscope %s: %.1f°\n',...
                 xyzStr(stg.rot(r)), x.out(i,r));    
        end       
        fprintf('   -> Deviation from ideal alignment in %s: %0.1f °\n',...
                 xyzStr(crys.alignAx(1)),eps.opt(i,1)); 
         if size(eps.opt,2) == 2
             fprintf('   -> Deviation from ideal alignment in %s: %0.1f °\n',...
                 xyzStr(crys.alignAx(2)),eps.opt(i,2)); 
         end   
         fprintf('-------------------------------------------------------------\n');    
   case 'FIB_FEIHelios'       
        %Variable input
        x = varargin{1}; 
        i = varargin{2};
        FIB = varargin{3};
        %Calculate upper and lower trench length (y)
        FIB.ang(i).dR = x.out(i,FIB.axs.rot);
        FIB.ang(i).t0 = x.out(i,FIB.axs.tilt);
        FIB.ang(i).t52 = x.out(i,FIB.axs.tilt)+52;
        FIB.ang(i).t52inv = -x.out(i,FIB.axs.tilt)+52;
        if FIB.mode
               FIB.y(i).upper = cosd(FIB.ang(i).t0)*FIB.trench.z*sind(FIB.trench.ang)/...
                             sind(90-FIB.ang(i).t0-FIB.trench.ang); 
               FIB.y(i).lower = cosd(-FIB.ang(i).t0)*FIB.trench.z*sind(FIB.trench.ang)/...
                             sind(90+FIB.ang(i).t0-FIB.trench.ang);
            %Screen output
            fprintf('*** FIB - FEI Helios ***\n');
            fprintf('Apply:\n');   
            fprintf('   -> Relative rotation of %.1f°\n',FIB.ang(i).dR);
            fprintf('   -> Tilt at lift-out position: %.1f°\n',FIB.ang(i).t0);
            fprintf('   -> Tilt at trenching position: %.1f°\n',FIB.ang(i).t52);
            if FIB.ang(i).t52 > FIB.ang(i).t52inv + 10 %Suggest 180 deg rotation
                fprintf('   -> Alternative: Tilt at trenching position: %.1f° + 180° relative rotation\n',FIB.ang(i).t52inv);
            end
            fprintf('Trench lengths for %.1f µm trench depth (z) and %.1f° trench angle:\n',FIB.trench.z,FIB.trench.ang);
            fprintf('   -> Trench length (y) at ''downhill position'': %.1f µm\n',FIB.y(i).lower);
            fprintf('   -> Trench length (y) at ''uphill position'': %.1f µm\n',FIB.y(i).upper);
            fprintf('-------------------------------------------------------------\n');
        end
end
end
%% xyzStr - Return string of major coordinate axes
function str = xyzStr(vec)
%function str = xyzStr(vec)
%vec: vector3d [xvector | yvector | zvector]
    if vec.x == 1 && vec.y == 0 && vec.z == 0
       str = 'X-axis';                                                     %Assign axis string
    elseif vec.x == 0 && vec.y == 1 && vec.z == 0
       str = 'Y-axis';                                                     %Assign axis string
    elseif vec.x == 0 && vec.y == 0 && vec.z == 1
       str = 'Z-axis';                                                     %Assign axis string
    else
       str = num2str(vec.xyz);     
       str = str(str ~= ' ');
    end    
end
%% tileFigs - Tile figures accross screen
function tileFigs() 
%function tileFigs() 
%Tile all figures evenly spread accros the screen
%% Initialization
mon = 1;                                                                   %Choose monitor number
offset.l = 70; offset.r = 0; offset.b = 70; offset.t = 0;                  %Offsets left right botton top (possible taskbars)
grid = [2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4; 
3 3 3 3 3 3 3 3 4 4 4 5 5 5 5 5 5 5 5 6]';                                 %Define figure grid
%% Find figures and screen dimension
h.figs = flip(findobj('type','figure'));                                   %Get figure handles
set(h.figs,'unit','pixels');                                               %Set figure units to [pxs]
nFigs = size(h.figs,1);                                                    %Get number of figures

scr.Sz = get(0,'MonitorPositions');                                        %Get screen size
scr.h = scr.Sz(mon,4)-offset.t;                                            %Get screen height
scr.w = scr.Sz(mon,3)-offset.l-offset.r;                                   %Get screen width
scr.orX = scr.Sz(mon,1)+offset.l;                                          %Get screen origin X
scr.orY = scr.Sz(mon,2);                                                   %Get screen origin Y
%% Check limits
if ~nFigs; warning('figures are not found'); return; end                   %Stop for no figures
if nFigs > 20; warning('too many figures(maximum = 20)'); return; end      %Check for limit of 20 figures
%% Define grid according to screen aspect ratio
if scr.w > scr.h %Widescreen
    n.h = grid(nFigs,1);                                                   %Define number of figures in height                                                 
    n.w = grid(nFigs,2);                                                   %Define number of figures in width         
else 
    n.h = grid(nFigs,2);                                                   %Define number of figures in height 
    n.w = grid(nFigs,1);                                                   %Define number of figures in width  
end 
%% Determine height and width for each figure
fig.h = (scr.h-offset.b)/n.h;                                              %Figure height
fig.w =  scr.w/n.w;                                                        %Figure width
%% Resize figures
k = 1;                                                                     %Initialize figure counter
for i =1:n.h %Loop over height
    for j = 1:n.w  %Loop over width
        if k > nFigs; return; end                                          %Stop when all figures have been resized 
        fig_pos = [scr.orX + fig.w*(j-1) scr.h-fig.h*i fig.w fig.h];   %Compute new figure position 
        set(h.figs(k),'OuterPosition',fig_pos);                            %Set new figure position
        k = k + 1;                                                         %Increase figure counter
    end 
end
end

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
