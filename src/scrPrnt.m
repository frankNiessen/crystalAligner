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
                xyzStr(crys.alignAx(i)),crys.cs{i},crys.type{i},reshape(char(crys.strMil{i})',1,[]),sym);
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
        i = varargin{3};
        %Screen output
        fprintf('\n*************************************************************');
        fprintf('\n Optimization problem %.0f - Parallel alignment of:\n',i);
        fprintf('   -> %s %s %s with microscope %s\n',crys.cs{1},crys.strMil{1}{i},crys.type{1},xyzStr(crys.alignAx(1)));
        if strcmp(optim.Alg,'gamultiobj')
           if size(crys.Miller{2},1)==1
               fprintf('   -> %s %s %s with microscope %s\n',crys.cs{2},char(crys.strMil{2}),crys.type{2},xyzStr(crys.alignAx(2)));
           else
               fprintf('   -> Either of %s %ss with microscope %s\n',crys.oMil{2}{1}.CS.LaueName,crys.type{2},char(crys.strMil{2}),xyzStr(crys.alignAx(2)));
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
                 xyzStr(stg.rot(r)), x.out(i,r)*stg.sign(r));
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
            fprintf('   -> Trench length (y) at ''downhill position'': %.1f °m\n',FIB.y(i).lower);
            fprintf('   -> Trench length (y) at ''uphill position'': %.1f °m\n',FIB.y(i).upper);
            fprintf('-------------------------------------------------------------\n');
        end
end
end