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
        for nrMill = 1:optim.order
           if crys{nrMill}.Miller.opt.useSym
               sym='(Symmetry activated)'; 
           else
               sym = '(Symmetry deactivated)'; 
           end
           fprintf(' - Microscope %s with "%s" %s %s\n',...
                crys{nrMill}.alignAx.char,crys{nrMill}.CS.mineral,...
                crys{nrMill}.Miller.char,sym);
        end
        fprintf('Rotational degrees of freedom:\n');
        for r = 1:size(stg.rot,1)
            fprintf(' - Rotation around microscope %s: %.1f to %.1f°\n',...
                    xyzStr(stg.rot(r)),stg.LB(r),stg.UB(r));
        end
        fprintf('-------------------------------------------------------------');
    case 'Optim'
        %Variable Input
        crys = varargin{1};
        optim = varargin{2};
        nrMill = varargin{3};
        %Screen output

        fprintf('\n*************************************************************\n');
        fprintf('Optimization problem %d - Parallel alignment of:\n',nrMill);
        fprintf('   -> %s %s with microscope %s\n',crys{1}.CS.mineral,crys{1}.Miller(nrMill).char, ...
                                                   crys{1}.alignAx.char);
        if strcmp(optim.Alg,'gamultiobj')
           if length(crys{2}.Miller)==1
               fprintf('   -> %s %s with microscope %s\n',crys{2}.CS.mineral,crys{2}.Miller.char,...
                                                          crys{2}.alignAx.char);
           else
               fprintf('   -> Either of %s %s with microscope %s\n',crys{2}.CS.mineral,crys{2}.Miller.char,...
                                                                    crys{2}.alignAx.char);
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
            fprintf('under consideration of weighting factors %.0f and %.0f for alignment with microscope axes %s and %s\n',...
                     optim.wghtFac(1),optim.wghtFac(2),crys{1}.alignAx.char,crys{2}.alignAx.char);
        end
   case 'Result'
        %Variable input
        dir = varargin{1};
        stg = varargin{2};
        x = varargin{3};
        eps = varargin{4};
        nrMill = varargin{5};
        crys = varargin{6};
        optim = varargin{7};
        %Screen Output
        fprintf('\n-------------------------------------------------------------\n');
        fprintf('*** Optimization results ***\n');
        if strcmp(optim.Alg,'ga')
            fprintf('Optimal parallel alignment of\n\t-Microscope axis %s with "%s" %s:\n',...
                     crys{1}.alignAx.char,crys{1}.CS.mineral,crys{1}.Miller(nrMill).char);
        elseif strcmp(optim.Alg,'gamultiobj')
            fprintf('\nOptimal parallel alignment of\n\t- microscope axis %s with "%s" %s and \n\t- microscope axis %s with "%s" %s:\n',...
                     crys{1}.alignAx.char,crys{1}.CS.mineral,crys{1}.Miller(nrMill).char,...
                     crys{2}.alignAx.char,crys{2}.CS.mineral,crys{2}.Miller.char);
        end

        for r = 1:size(stg.rot,1)
            fprintf('   -> Rotation around microscope %s: %.2f°\n',...
                 xyzStr(stg.rot(r)), x.out(nrMill,r));
        end
        fprintf('   -> Deviation from ideal alignment in %s: %0.2f°\n',...
                 crys{1}.alignAx.char,eps.opt(nrMill,1));
         if size(eps.opt,2) == 2
             fprintf('   -> Deviation from ideal alignment in %s: %0.2f°\n',...
                 crys{2}.alignAx.char,eps.opt(nrMill,2));
         end
         fprintf('-------------------------------------------------------------\n');
   case 'FIB_FEIHelios'
        %Variable input
        x = varargin{1};
        nrMill = varargin{2};
        FIB = varargin{3};
        %Calculate upper and lower trench length (y)
        FIB.ang(nrMill).dR = x.out(nrMill,FIB.axs.rot);
        FIB.ang(nrMill).t0 = x.out(nrMill,FIB.axs.tilt);
        FIB.ang(nrMill).t52 = x.out(nrMill,FIB.axs.tilt)+52;
        FIB.ang(nrMill).t52inv = -x.out(nrMill,FIB.axs.tilt)+52;
        if FIB.mode
               FIB.y(nrMill).upper = cosd(FIB.ang(nrMill).t0)*FIB.trench.z*sind(FIB.trench.ang)/...
                             sind(90-FIB.ang(nrMill).t0-FIB.trench.ang);
               FIB.y(nrMill).lower = cosd(-FIB.ang(nrMill).t0)*FIB.trench.z*sind(FIB.trench.ang)/...
                             sind(90+FIB.ang(nrMill).t0-FIB.trench.ang);
            %Screen output
            fprintf('*** FIB - FEI Helios ***\n');
            fprintf('Apply:\n');
            fprintf('   -> Relative rotation of %.2f°\n',FIB.ang(nrMill).dR);
            fprintf('   -> Tilt at lift-out position: %.2f°\n',FIB.ang(nrMill).t0);
            fprintf('   -> Tilt at trenching position: %.2f°\n',FIB.ang(nrMill).t52);
            if FIB.ang(nrMill).t52 > FIB.ang(nrMill).t52inv + 10 %Suggest 180 deg rotation
                fprintf('   -> Alternative: Tilt at trenching position: %.2f° + 180° relative rotation\n',FIB.ang(nrMill).t52inv);
            end
            fprintf('Trench lengths for %.1f µm trench depth (z) and %.2f° trench angle:\n',FIB.trench.z,FIB.trench.ang);
            fprintf('   -> Trench length (y) at ''downhill position'': %.1f µm\n',FIB.y(nrMill).lower);
            fprintf('   -> Trench length (y) at ''uphill position'': %.1f µm\n',FIB.y(nrMill).upper);
            fprintf('-------------------------------------------------------------\n');
        end
end
end