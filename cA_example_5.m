% *************************************************************************
% crystalAligner.m - Crystal Alignment Tool for Electron Microscopy
% .........................................................................
% Example 5
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
% 2021/04, Improvements implemented by Dr. Vivian Tong, NPL
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
clear variables; close all; home; 
fprintf('\n*************************************************************');
fprintf('\n               crystalAligner Example 5\n');
fprintf('*************************************************************\n\n');
fprintf('-> Starting up MTEX ...\n');                                      %ScreenPrint
currentFolder;                                                             %Add subfolder to path
iniExtLibs;                                                                %Initialize MTEX and check for MATLAB toolboxes

%% Definition of alignment problem
% ****************************** Crystals *********************************
% *** Crystal Alignment Objective 1
crys(1).SS = specimenSymmetry('triclinic');                                %Define specimen symmetry ['triclinic'|'orthorhombic']
crys(1).CS = crystalSymmetry('m-3m',[3.24,3.24,3.24],'mineral','beta');   %Define crystal symmetry (https://mtex-toolbox.github.io/CrystalSymmetries.html)
crys(1).ori = orientation('Euler',246*degree,36*degree,75*degree, ...
                                  crys(1).CS,crys(1).SS);                  %Define crystal orientation (https://mtex-toolbox.github.io/OrientationDefinition.html)   
crys(1).alignAx = xvector;                                                 %Microscope axis for alignment with crystal direction/plane; Examples: zvector; [.5 .5 1]; xvector; ...
crys(1).Miller = Miller(1,1,0,crys(1).CS,crys(1).SS,'hkl');                %Define at least one crystal direction (https://mtex-toolbox.github.io/CrystalDirections.html)
crys(1).Miller.opt.useSym = true;                                          %Apply crystal symmetry: true/false
% *** Crystal Alignment Objective 2
crys(2).SS = specimenSymmetry('triclinic');                                %Define specimen symmetry ['triclinic'|'orthorhombic']
crys(2).CS = crystalSymmetry('P63/mmc',[2.95,2.95,4.686],'mineral','alphaP'); %Define crystal symmetry (https://mtex-toolbox.github.io/CrystalSymmetries.html)
crys(2).ori = orientation('Euler',91*degree,94*degree,13*degree, ...
                                  crys(2).CS,crys(2).SS);                  %Define crystal orientation (https://mtex-toolbox.github.io/OrientationDefinition.html)       
crys(2).alignAx = zvector;                                                 %Microscope axis for alignment with crystal direction/plane; Examples: zvector; [.5 .5 1]; xvector; ...
crys(2).Miller  = Miller(1,1,-2,0,crys(2).CS,crys(2).SS,'uvtw');           %Define one crystal direction (https://mtex-toolbox.github.io/CrystalDirections.html)
crys(2).Miller.opt.useSym = true;                                          %Apply crystal symmetry: true/false
% ******************************* Stage ***********************************                                          
stg.rot     = [xvector; -zvector];                                          %Stage rotation axes                                        
stg.LB      = [    0     -180  ];                                          %Lower bound [°]
stg.UB      = [   20      180  ];                                          %Upper bound [°]                                         
stg.order   = [    2      1    ];                                          %Hierarchy / order of rotation: Rotation 1 before 2 before 3; Example: [3 1 2];
%Genetic algorithm
optim.order = length(crys);                                                %Order of optimization problem (do not change)
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

%% Setup Optimization options
crys = checkerror(crys);                                                   %Error checking
optim = setOptimOpts(optim);                                               %Optimization initialization function
scrPrnt('Ini',crys,stg,optim);                                             %Screen print of optimization objectives and limits parameters

%% Run crystalAligner
[oNew,stgRot,x,eps] = crystalAligner(crys,stg,optim,FIB);
fprintf('\n -> All done!\n');                                       

