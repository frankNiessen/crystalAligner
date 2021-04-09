function stereoProj(ori,Miller,titleStr)
%function stereoProj(o,Dir,Dirstrs,titleStr)
%Plot stereographic projection of all equivalent crystal directions to
%Miller crystal direction 'Dir' with respect crystal orientation 'o'.
%'Dirstrs' are legend strings and 'titleStr' defines the title string
%% Initialization
markerSz = 12;                                                             %Define marker size
markers = {'o','d','s','v','h','^'};                                       %Define marker style
colors = {'k','r','b','g','m','c'};                                        %Define marker colors
lgFntSz = 20;                                                              %FontSize for legend
%% Plotting
figure;                                                                    %Create new figure
for ii = 1:length(Miller) %Loop over crystal directions
    if Miller.opt.useSym
        dirs = symmetrise(Miller(ii));
    else
        dirs = Miller(ii);
    end
    plot(ori*dirs,'antipodal','grid','backgroundcolor','w',...
                 'MarkerSize',markerSz,'MarkerEdgeColor','k',...
                 'MarkerFaceColor',colors{ii},'Marker',markers{ii});         %Plot equivalent crystal directions
    hold on                                                                %Allow additive plotting
    annotate([xvector,yvector],'label',{'X','Y'},'backgroundcolor','w');   %Annotate labels
    L(ii) = plot(nan(6,1),nan(6,1),'MarkerSize',markerSz,...
           'MarkerEdgeColor','k','MarkerFaceColor',colors{ii},...
           'Marker',markers{ii},'Color','w');                               %Create 'Ghost plot' for creating a nice legend
    legStr{ii} = Miller(ii).char;
end
%% Finishing
l = legend(L,legStr,'FontName','TimesNewRoman','FontSize',20);            %Create legend
l.FontSize = lgFntSz;                                                      %Change legend fontsize
set(gcf,'name',titleStr);
end