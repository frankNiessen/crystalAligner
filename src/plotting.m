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