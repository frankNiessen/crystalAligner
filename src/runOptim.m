%% runOptim - Optimization function
function [oNew,stgRot,x,eps] = runOptim(crys,stg,optim,FIB)
ori = {cell2mat(crys).ori};
alignAx = [cell2mat(crys).alignAx];
for ii = 1:length(crys{1}.Miller) %Loop over crystal directions
    scrPrnt('Optim',crys,optim,ii);                                     %Screen output optimization problem
    % *** Multiobjective opimization **************************************
    if crys{1}.Miller.opt.useSym
        dirs = symmetrise(crys{1}.Miller(ii));
    else
        dirs = crys{1}.Miller(ii);
    end 

    if strcmp(optim.Alg,'gamultiobj')
        % *** Define optimization function *****

        fMin = @(x)minFunc2D(x,ori,dirs,crys{2}.Miller,alignAx,stg);
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
        x.opt(ii,:) = x.ga(ind,:);                                         %Save optimal solution
        eps.opt(ii,:) = eps.ga(ind,:);                                     %Save misalignment of optimal solution
        x.out(ii,:) = x.opt(ii,:).*stg.sign;                               %Adapt possible different stage rotation convention for output

        %Plot optimal solution
        if optim.plot
            h.ax = findobj('type','axes','tag','gaplotpareto');            %Find axes
            hold(h.ax,'on');                                               %Hold on
            plot([0,eps.opt(ii,1)],[0,eps.opt(ii,2)],'k');                   %Plot line
            plot(eps.opt(ii,1),eps.opt(ii,2),'*k');                          %Plot marker
        end
        %Find out which second crystal direction was aligned
        rotTot = rotation('Euler',0,0,0);
        for r = 1:size(stg.order,2)
            axNr = stg.order(r);
            rot(axNr) = rotation('axis',stg.rot(axNr),'angle',x.opt(axNr)*degree);     %Compute rotation around microscope rotation axis 'r' 'stg.rot(r)'
            rotTot = rot(axNr)*rotTot;
        end
        for jj=1:length(crys{2}.Miller)
            if crys{2}.Miller.opt.useSym
                dirs = symmetrise(crys{2}.Miller(jj));
            else
                dirs = crys{2}.Miller(jj);
            end
           epsTemp(jj)=min(angle(rotTot*ori{2}*dirs,alignAx(2))/degree); %Find misalignment of microscope axis 2 'crys.alignAx(2)' with 'o*CD' subject to rotations 'rotX*rotZ'
        end
        [~,crystDirInd] = min(epsTemp);
        crys{2}.optAx = crys{2}.Miller(crystDirInd);
    % *** Singleobjective opimization *************************************
    elseif strcmp(optim.Alg,'ga')
        % *** Define optimization function *****
        fMin = @(x)minFunc(x,ori{1},dirs,stg,crys{1}.alignAx);
        % *** Run optimization *****
        [x.opt(ii,:),eps.opt(ii,:)] = ga(fMin,size(stg.rot,1),[],[],[],[],stg.LB,stg.UB,[],optim.opt); %genetic algorithm
         x.out(ii,:) = x.opt(ii,:).*stg.sign;                                %Adapt possible different stage rotation convention for output

    else
        error('Invalid choice of optimization algorithm');                 %No valid choice of alcrys.oithm
    end
                           %Optimized crystal direction with second axis
    %Output
    scrPrnt('Result',dir,stg,x,eps,ii,crys,optim);                          %General results screen output
    scrPrnt('FIB_FEIHelios',x,ii,FIB);                                      %Instrument specific screen output
    rotTot = rotation('Euler',0,0,0);                                                            %Ini
    for r = 1:size(stg.rot,1)
        axNr = stg.order(r);
        stgRot{ii}.ax{axNr} = rotation('axis',stg.rot(axNr),...
                                       'angle',x.opt(ii,axNr)*degree);     %Convert rotation around microscope axis r 'stg.rot(r)' to Euler angles
        rotTot = stgRot{ii}.ax{axNr}*rotTot;                               %Total rotation
    end
    for c = 1:length(ori)
        oNew{ii}{c} = rotTot*ori{c};                                       %Compute new crystal orientation after applied stage tilt and rotation
    end
end
end