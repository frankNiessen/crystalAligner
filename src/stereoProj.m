function stereoProj(o,Dir,Dirstrs,sym,titleStr)
%function stereoProj(o,Dir,Dirstrs,titleStr)
%Plot stereographic projection of all equivalent crystal directions to
%Miller crystal direction 'Dir' with respect crystal orientation 'o'.
%'Dirstrs' are legend strings and 'titleStr' defines the title string
%% Initialization
figure;                                                                    %Create new figure
markerSz = 12;                                                             %Define marker size
markers = {'o','d','s','v','h','^'};                                       %Define marker style
colors = {'k','r','b','g','m','c'};                                        %Define marker colors
lgFntSz = 20;                                                              %FontSize for legend
%% Plotting
for i = 1:length(Dir) %Loop over crystal directions
    if sym
        dirs = symmetrise(Dir{i});
    else
        dirs = Dir{i};
    end
    plot(o*dirs,'antipodal','grid','backgroundcolor','w',...
                 'MarkerSize',markerSz,'MarkerEdgeColor','k',...
                 'MarkerFaceColor',colors{i},'Marker',markers{i});         %Plot equivalent crystal directions
    hold on                                                                %Allow additive plotting
    annotate([xvector,yvector],'label',{'X','Y'},'backgroundcolor','w');   %Annotate labels
    L(i) = plot(nan(6,1),nan(6,1),'MarkerSize',markerSz,...
           'MarkerEdgeColor','k','MarkerFaceColor',colors{i},...
           'Marker',markers{i},'Color','w');                               %Create 'Ghost plot' for creating a nice legend
end
%% Finishing
l = legend(L,Dirstrs,'FontName','TimesNewRoman','FontSize',20);            %Create legend
l.FontSize = lgFntSz;                                                      %Change legend fontsize
set(gcf,'name',titleStr);
end