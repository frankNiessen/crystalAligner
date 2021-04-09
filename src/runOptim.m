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
                    'Choose a Pareto-solution: ax{1}[ï¿½] ax{2}[ï¿½] eps1[ï¿½] eps2[ï¿½]:',...
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