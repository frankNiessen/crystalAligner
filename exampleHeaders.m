%% Example E1
% *** Crystal Alignment Objective 1
crys.o(1,:)     = [313 15 137];                                            %Crystal orientation in Euler angles [pih1 Phi phi2]       
crys.cs{1}      = 'cubic';                                                 %Crystal structure string (follow MTEX convention)
crys.alignAx(1) = zvector;                                                 %Microscope axis for alignment with crystal direction/plane; Examples: zvector; [.5 .5 1]; xvector; ...
crys.Miller{1}  = [1 1 1; 1 0 0];                                          %Miller indices for alignment (in Multiobjective Optimization several Miller-sets will start several optimizations);  Examples: [1 0 0; 1 1 0]; [-1 2 1]; ...
crys.type{1}    = 'hkl';                                                   %Type of Miller: 'hkl': Crystal plane; 'uvw': Crystal direction
crys.sym(1)     = 1;                                                       %Apply crystal symmetry: 1: yes 0: no
% ******************************* Stage ***********************************                                          
stg.rot     = [xvector; zvector];                                          %Stage rotation axes                                        
stg.LB      = [    0     -180  ];                                          %Lower bound [°]
stg.UB      = [   20      180  ];                                          %Upper bound [°]                                         
stg.sign    = [    1      -1   ];                                          %Sign 1: Right hand rule convention; -1: Left hand rule convention
stg.order   = [    2      1    ];                                          %Hierarchy / order of rotation: Rotation 1 before 2 before 3; Example: [3 1 2];
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
FIB.mode = 0;                                                              %Flag: FIB liftout output [1|0]
FIB.trench.ang = 52;                                                       %Trenching - or look-in - angle of Trench [°]
FIB.trench.z = 15;                                                         %Trench depth 'z' [µm]
FIB.axs.tilt = 1;                                                          %Index of tilt axis in 'axs.rot'
FIB.axs.rot = 2;                                                           %Index of rotation axis in 'axs.rot'                                                        
%Output
optim.plot = 1;                                                            %Plotting 1: On 0: Off

%% Example E2
% ****************************** Crystals *********************************
% *** Crystal Alignment Objective 1
crys.o(1,:)     = [61 42 9];                                               %Crystal orientation in Euler angles [pih1 Phi phi2]       
crys.cs{1}      = 'cubic';                                                 %Crystal structure string (follow MTEX convention)
crys.alignAx(1) = yvector;                                                 %Microscope axis for alignment with crystal direction/plane; Examples: zvector; [.5 .5 1]; xvector; ...
crys.Miller{1}  = [1 1 3];                                                 %Miller indices for alignment (in Multiobjective Optimization several Miller-sets will start several optimizations);  Examples: [1 0 0; 1 1 0]; [-1 2 1]; ...
crys.type{1}    = 'hkl';                                                   %Type of Miller: 'hkl': Crystal plane; 'uvw': Crystal direction
crys.sym(1)     = 1;                                                       %Apply crystal symmetry: 1: yes 0: no
% *** Crystal Alignment Objective 2
crys.o(2,:)     = [61 42 9];                                               %Crystal orientation in Euler angles [pih1 Phi phi2]       
crys.cs{2}      = 'cubic';                                                 %Crystal structure string (follow MTEX convention)
crys.alignAx(2) = zvector;                                                 %Microscope axis for alignment with crystal direction/plane; Examples: zvector; [.5 .5 1]; xvector; ...
crys.Miller{2}  = [1 1 0];                                                 %Miller indices for alignment (in Multiobjective Optimization several Miller-sets will start several optimizations);  Examples: [1 0 0; 1 1 0]; [-1 2 1]; ...
crys.type{2}    = 'hkl';                                                   %Type of Miller: 'hkl': Crystal plane; 'uvw': Crystal direction
crys.sym(2)     = 1;                                                       %Apply crystal symmetry: 1: yes 0: no
% ******************************* Stage ***********************************                                          
stg.rot     = [xvector; zvector];                                          %Stage rotation axes                                        
stg.LB      = [    0     -180  ];                                          %Lower bound [°]
stg.UB      = [   20      180  ];                                          %Upper bound [°]                                         
stg.sign    = [    1      -1   ];                                          %Sign 1: Right hand rule convention; -1: Left hand rule convention
stg.order   = [    2      1    ];                                          %Hierarchy / order of rotation: Rotation 1 before 2 before 3; Example: [3 1 2];
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

%% Example E3
% ****************************** Crystals *********************************
% *** Crystal Alignment Objective 1
crys.o(1,:)     = [261 43 28];                                             %Crystal orientation in Euler angles [pih1 Phi phi2]       
crys.cs{1}      = 'cubic';                                                 %Crystal structure string (follow MTEX convention)
crys.alignAx(1) = xvector;                                                 %Microscope axis for alignment with crystal direction/plane; Examples: zvector; [.5 .5 1]; xvector; ...
crys.Miller{1}  = [0 1 1];                                                 %Miller indices for alignment (in Multiobjective Optimization several Miller-sets will start several optimizations);  Examples: [1 0 0; 1 1 0]; [-1 2 1]; ...
crys.type{1}    = 'hkl';                                                   %Type of Miller: 'hkl': Crystal plane; 'uvw': Crystal direction
crys.sym(1)     = 1;                                                       %Apply crystal symmetry: 1: yes 0: no
% *** Crystal Alignment Objective 2
crys.o(2,:)     = [175 20 102];                                            %Crystal orientation in Euler angles [pih1 Phi phi2]       
crys.cs{2}      = 'orthorhombic';                                          %Crystal structure string (follow MTEX convention)
crys.alignAx(2) = zvector;                                                 %Microscope axis for alignment with crystal direction/plane; Examples: zvector; [.5 .5 1]; xvector; ...
crys.Miller{2}  = [0 0 1];                                                 %Miller indices for alignment (in Multiobjective Optimization several Miller-sets will start several optimizations);  Examples: [1 0 0; 1 1 0]; [-1 2 1]; ...
crys.type{2}    = 'hkl';                                                   %Type of Miller: 'hkl': Crystal plane; 'uvw': Crystal direction
crys.sym(2)     = 0;                                                       %Apply crystal symmetry: 1: yes 0: no
% ******************************* Stage ***********************************                                          
stg.rot     = [xvector; zvector];                                          %Stage rotation axes                                        
stg.LB      = [    0     -180  ];                                          %Lower bound [°]
stg.UB      = [   20      180  ];                                          %Upper bound [°]                                         
stg.sign    = [    1      -1   ];                                          %Sign 1: Right hand rule convention; -1: Left hand rule convention
stg.order   = [    2      1    ];                                          %Hierarchy / order of rotation: Rotation 1 before 2 before 3; Example: [3 1 2];
% ************************* Genetic algorithm *****************************
%Genetic algorithm
optim.popSz = 500;                                                         %Population size
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
FIB.mode = 0;                                                              %Flag: FIB liftout output [1|0]
FIB.trench.ang = 52;                                                       %Trenching - or look-in - angle of Trench [°]
FIB.trench.z = 15;                                                         %Trench depth 'z' [µm]
FIB.axs.tilt = 1;                                                          %Index of tilt axis in 'axs.rot'
FIB.axs.rot = 2;                                                           %Index of rotation axis in 'axs.rot'                                                        
%Output
optim.plot = 1;                                                            %Plotting 1: On 0: Off

%% Example E4
% *** Crystal Alignment Objective 1
crys.o(1,:)     = [313 15 137];                                            %Crystal orientation in Euler angles [pih1 Phi phi2]       
crys.cs{1}      = 'cubic';                                                 %Crystal structure string (follow MTEX convention)
crys.alignAx(1) = zvector;                                                 %Microscope axis for alignment with crystal direction/plane; Examples: zvector; [.5 .5 1]; xvector; ...
crys.Miller{1}  = [1 1 1; 1 0 0];                                          %Miller indices for alignment (in Multiobjective Optimization several Miller-sets will start several optimizations);  Examples: [1 0 0; 1 1 0]; [-1 2 1]; ...
crys.type{1}    = 'hkl';                                                   %Type of Miller: 'hkl': Crystal plane; 'uvw': Crystal direction
crys.sym(1)     = 1;                                                       %Apply crystal symmetry: 1: yes 0: no
% ******************************* Stage ***********************************                                          
stg.rot     = [xvector; yvector; zvector];                                 %Stage rotation axes                                        
stg.LB      = [    0     -45      -180  ];                                 %Lower bound [°]
stg.UB      = [   20      45       180  ];                                 %Upper bound [°]                                         
stg.sign    = [    1      1        -1   ];                                 %Sign 1: Right hand rule convention; -1: Left hand rule convention
stg.order   = [    3      1         2   ];                                 %Hierarchy / order of rotation: Rotation 1 before 2 before 3; Example: [3 1 2];
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
FIB.mode = 0;                                                              %Flag: FIB liftout output [1|0]
FIB.trench.ang = 52;                                                       %Trenching - or look-in - angle of Trench [°]
FIB.trench.z = 15;                                                         %Trench depth 'z' [µm]
FIB.axs.tilt = 1;                                                          %Index of tilt axis in 'axs.rot'
FIB.axs.rot = 2;                                                           %Index of rotation axis in 'axs.rot'                                                        
%Output
optim.plot = 1;                                                            %Plotting 1: On 0: Off