% *************************************************************************
% crystalAligner.m - Crystal Alignment Tool for Electron Microscopy
% .........................................................................
% Example 5b
% Concurrent alignment of the {011} plane normals and <11-20> crystal 
% directions in a cubic and a hexagonal crystal with the x-Axis and z-Axis 
% of the microscope, respectively.
% *************************************************************************
% If you use this script for research and find it useful I would appreciate a
% citation of the associated article:
% [crystalAligner: a computer program to align crystal directions in a scanning
% electron microscope by global optimization, F. Niessen, Journal of Applied
% Crystallography, 2020, doi: 10.1107/S1600576719017345]
% *************************************************************************
% Dr. Frank Niessen, University of Wollongong, Australia, 2018
% contactnospam@fniessen.com (remove the nospam to make this email address
% work)
% 2021/04/08, Improvements implemented by Dr. Vivian Tong, NPL
% *************************************************************************
% Execution of the software requires installation of
% * mtex (http://mtex-toolbox.github.io/)
% * MATLAB Optimization toolbox
% * MATLAB Global optimization toolbox
% * Optional: MATLAB Parallel Computing Toolbox
% Compatibility tested with MATLAB R2104a - R2016b and MTEX 4.2.1 to 5.6
% *************************************************************************
% Copyright 2018 Frank Niessen (see attached MIT license file)

%% Initialization
clear variables; close all; home; ; close all; home; 
fprintf('\n*************************************************************');
fprintf('\n                  crystalAligner v. 1.1 \n');
fprintf('*************************************************************\n\n');
fprintf('\n -> Starting up MTEX ...');                                     %ScreenPrint
currentFolder;                                                             %Add subfolder to path
iniExtLibs;                                                                %Initialize MTEX and check for MATLAB toolboxes

%% Definition of alignment problem
% ****************************** Crystals *********************************
crys.ss = 'orthorhombic';
% *** Crystal Alignment Objective 1
crys.o{1}     = [246 36 75];                                               %Crystal orientation in Euler angles [pih1 Phi phi2]       
crys.cs{1}      = crystalSymmetry('6/mmm', [2.9065 2.9065 2.8366], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'W C');                                                  %Crystal structure string (follow MTEX convention)
crys.alignAx(1) = xvector;                                                 %Microscope axis for alignment with crystal direction/plane; Examples: zvector; [.5 .5 1]; xvector; ...
crys.Miller{1}  = [1 0 -1 1; 1 1 -2 0];                                    %Miller indices for alignment (in Multiobjective Optimization several Miller-sets will start several optimizations);  Examples: [1 0 0; 1 1 0]; [-1 2 1]; ...
crys.type{1}    = 'hkil';                                                  %Type of Miller: 'hkl/hkil': Crystal plane; 'uvw/UVTW': Crystal direction
crys.sym(1)     = true;                                                    %Apply crystal symmetry: true/false
% *** Crystal Alignment Objective 2
crys.cs{2}      = 'cubic';                                               %Crystal structure string (follow MTEX convention)
crys.alignAx(2) = zvector;                                                 %Microscope axis for alignment with crystal direction/plane; Examples: zvector; [.5 .5 1]; xvector; ...
crys.Miller{2}  = [1 1 0];                                              %Miller indices for alignment (in Multiobjective Optimization several Miller-sets will start several optimizations);  Examples: [1 0 0; 1 1 0]; [-1 2 1]; ...
crys.o{2}       = orientation('Euler',20*degree,30*degree,10*degree,crys.cs{2},crys.ss);                                              %Crystal orientation in Euler angles [pih1 Phi phi2]       
crys.type{2}    = 'uvw';                                                  %Type of Miller: 'hkl': Crystal plane; 'uvw': Crystal direction
crys.sym(2)     = true;                                                    %Apply crystal symmetry: true/false
% ******************************* Stage ***********************************                                          
stg.rot     = [xvector; zvector];                                          %Stage rotation axes                                        
stg.LB      = [    0    57  ];                                             %Lower bound [°]
stg.UB      = [   20    57  ];                                             %Upper bound [°]                                         
stg.sign    = [    1      -1   ];                                          %Sign 1: Right hand rule convention; -1: Left hand rule convention
stg.order   = [    2      1    ];                                          %Hierarchy / order of rotation: Rotation 1 before 2 before 3; Example: [3 1 2];
% ************************* Genetic algorithm *****************************
%Genetic algorithm
optim.popSz = 200;                                                         %Population size
optim.funcTol = 0.1;                                                       %FunctionTolerance
optim.maxStallGen = 10;                                                    %Maximum stall generations
optim.iterOut = false;                                                     %Writing output for each iteration in subFolder 'iterOut': true/false
%Multiobjective genetic algorithm settings
optim.wghtFac = [1,1];                                                     %Weighting factors for TOPSIS multiobjective decision making method
optim.multiCore = false;                                                   %Utilization of parallel processing (switch off if errors ocur): true/false
optim.hybridFcn = false;                                                   %Flag: Use a hybrid function to (may speed op convergence but compromise diversity of solution space): true/false
optim.autoSol = true;                                                      %Flag: Pick optimum solution automatically by distance of Pareto solution from the optimal solution: true/false
% ***************************** Optional **********************************
%FIB liftout calculations
FIB.mode = false;                                                          %Flag: FIB liftout output: true/false
FIB.trench.ang = 52;                                                       %Trenching - or look-in - angle of Trench [°]
FIB.trench.z = 15;                                                         %Trench depth 'z' [µm]
FIB.axs.tilt = 1;                                                          %Index of tilt axis in 'axs.rot'
FIB.axs.rot = 2;                                                           %Index of rotation axis in 'axs.rot'                                                        
%Output
optim.plot = true;                                                         %Plotting 1: true/false

%% Error checking and PreProcessing
crys = checkerror(crys);                                                   %Error checking
%PreProcessing
crys.nrObj = length(crys.Miller);                                          %Number of alignment objectives
if ~isfield(crys,'ss') %if it wasn't already defined
    crys.ss = 'orthorhombic';                                                  %Specimen system
end

%% Setup Optimization options
optim = setOptimOpts(crys,optim);                                          %Optimization initialization function

%% Define Crystal system, orientation and crystal directions in mtex & Plot orientations
[ori,CS,SS,crys] = defOri(crys,optim);                                     %Orientation definition function
scrPrnt('Ini',crys,stg,optim);                                             %Screen print of optimization objectives and limits parameters

%% Optimization - Multiobjective genetic algorithm
[oNew,stgRot,x,eps] = runOptim(crys,stg,optim,ori,FIB,SS);                 %Optimization function

%% Plot stereographic projection and inverse polefigure of aligned equivalent crystal directions
plotting(optim,oNew,crys);                                                 %Plotting function
epsilon = eps.opt;                                                         %Assign return value
rot = x.opt;                                                               %Assign return value
fprintf('\n -> All done!\n');                                              %ScreenPrint
