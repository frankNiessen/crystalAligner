function plotting(optim,oNew,crys)
if optim.plot
    for cc = 1:optim.order
        for ii = 1:length(crys{cc}.Miller)
            titleStr = sprintf('%s aligned with "%s" %s',crys{cc}.alignAx.char, crys{cc}.CS.mineral, crys{cc}.Miller(ii).char);
            stereoProj(oNew{ii}{cc},crys{cc}.Miller,titleStr);%Plot new crystal directions of interest for z-axis
            figure;                                                        %Create figure
            plotIPDF(oNew{ii}{cc},[xvector,yvector,zvector],'antipodal',...
                                  'MarkerFaceColor','k', 'figsize','small'); %Plot inverse pole-figure (IPF) [x:(100), y:(010) and z:(001)]
            set(gcf,'name',titleStr);                                      %set title
        end
    end
    tileFigs();                                                                %Align figures
end
end