function plotOrientations(optim,ori,crys,mode)
if isa(ori{1},'orientation'); ori = {ori}; end                             %Unifying datastructure

if optim.plot
    for ii = 1:length(ori)   
        for cc = 1:length(ori{ii})
            if strcmpi(mode,'result')
               if cc == 2 
                   titleStr = sprintf('OPTIM %d OBJECTIVE %d "%s": %s aligned with %s',...
                              ii,cc,crys(cc).CS.mineral,crys(cc).Miller.char,crys(cc).alignAx.char); 
               else
                   titleStr = sprintf('OPTIM %d OBJECTIVE %d "%s": %s aligned with %s',...
                              ii,cc,crys(cc).CS.mineral,crys(cc).Miller(ii).char,crys(cc).alignAx.char);
               end
            elseif strcmpi(mode,'initial')
                titleStr = sprintf('INITIAL "%s": %s to be aligned with %s',...
                           crys(cc).CS.mineral,crys(cc).Miller.char,crys(cc).alignAx.char);
            else
                return
            end
            figure;
            stereoProj(ori{ii}{cc},crys(cc).Miller,titleStr);          %Plot new crystal directions of interest for z-axis
            figure;                                                
            plotIPDF(ori{ii}{cc},[xvector,yvector,zvector],'antipodal',...
                              'MarkerFaceColor','k', 'figsize','small');%Plot inverse pole-figure (IPF) [x:(100), y:(010) and z:(001)]
            set(gcf,'name',titleStr);                                  %set title    
        end       
    end
    tileFigs();                                                            %Align figures
end
end

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