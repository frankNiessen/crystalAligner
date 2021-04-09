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
err = min(angle(rotTot*o{1}*CD,crys.alignAx(1))/degree);                   %Find misalignment of microscope axis 1 'crys.alignAx(1)' with 'o*CD' subject to rotation 'rotTot'
end