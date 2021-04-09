function eps = minFunc2D(x,ori,CD1,CD2,alignAx,stg)
%function err = minFunc2D(x,o,CD,tensD,wghtFac)
%Objective function calculating the misalignment of the microscope
%axes 1 and 2 in 'crys.alignAx' with crystal directions 'CD1' and 'CD2'
%rotated by tilt/rotation angles around microscope axes 1 and 2 in 'stg.rot'
%under consideration of weighting factors 'wFac'

rotTot = rotation('Euler',0,0,0);
for r = 1:size(stg.order,2)
    axNr = stg.order(r);
    rot(axNr) = rotation('axis',stg.rot(axNr),'angle',x(axNr)*degree);     %Compute rotation around microscope rotation axis 'r' 'stg.rot(r)'
    rotTot = rot(axNr)*rotTot;
end
eps(1) = min(angle(rotTot*ori{1}*CD1,alignAx(1))/degree);                  %Find misalignment of microscope axis 1 'crys.alignAx(1)' with 'o*CD' subject to rotation 'rotTot'
epsTemp = zeros(length(CD2),1);                                            %Preallocate memory for 'epsTemp'
for i=1:length(CD2)
    if CD2.opt.useSym
        dirs = symmetrise(CD2(i));
    else
        dirs = CD2(i);
    end
   epsTemp(i)=min(angle(rotTot*ori{2}*dirs,alignAx(2))/degree);            %Find misalignment of microscope axis 2 'crys.alignAx(2)' with 'o*CD' subject to rotations 'rotX*rotZ'
end
eps(2)=min(epsTemp);                                                       %Take the smallest misalignment in axis 2 as optimum solution
end