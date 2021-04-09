function [o,cs,ss,crys] = defOri(crys,optim)
fprintf('   -> Defining crystal system, orientation and directions ...\n');  %ScreenPrint

% for c = 1:optim.order %Loop over crystals
%     %Define Specimen system
%     if isa(crys.ss,'crystalSymmetry')
%         ss = crys.ss;
%         crys.ss=ss.lattice.char; %overwrite crys.cs to maintain compatibility
%     elseif isa(crys.ss,'char')
%         ss = crystalSymmetry(crys.ss);
%     else
%         warning(['Wrong input type for specimen symmetry in crystal '...
%             num2str(c) '? I should be either a character string or a mtex crystalSymmetry object.'...
%             'using default orthorhombic specimen symmetry!']);
%         ss = crystalSymmetry('orthorhombic');
%     end
%     
%     %Define Crystal system
%     if isa(crys.cs{c},'crystalSymmetry')
%         cs{c} = crys.cs{c};
%         crys.cs{c}=cs{c}.lattice.char; %overwrite crys.cs to maintain compatibility
%     elseif isa(crys.cs{c},'char')
%         cs{c} = crystalSymmetry(crys.cs{c});
%     else
%         warning(['Wrong input type for crystal symmetry in crystal '...
%             num2str(c) '? I should be either a character string or a mtex crystalSymmetry object.'...
%             'using default cubic crystal symmetry!']);
%         cs{c} = crystalSymmetry('cubic');
%     end
%     
%     %Define crystal orientation
%     if isa(crys.o,'cell') && isa(crys.o{c},'orientation')                                  % mtex orientation object - direct copy
%         o{c} = crys.o{c};
%         crys.o{c}=[o{c}.phi1 o{c}.Phi o{c}.phi2];                                           %overwrite crys.o to maintain compatibility
%     elseif isa(crys.o,'cell') && isa(crys.o{c},'numeric') && length(crys.o{c})==3          % cell array of euler angle triplets
%         o{c} =  orientation('Euler',crys.o{c}(1)*degree,crys.o{c}(2)*degree,...
%             crys.o{c}(3)*degree,cs{c},ss);
%     elseif isa(crys.o,'numeric')  && size(crys.o,2)==3                                     % numeric array of euler angles
%         o{c} =  orientation('Euler',crys.o(c,1)*degree,crys.o(c,2)*degree,...
%             crys.o(c,3)*degree,cs{c},ss);              %Compute crystal orientation
%     else
%         warning(['Wrong input type for orientation in crystal '...
%             num2str(c) '? I should be a cell array of 3-number vectors (eul angs), '...
%             'a c*3 matrix of numbers (eul angs), or a mtex orientation object'.'...
%             'using default euler angles [0 0 0]!']);
%         o{c} =  orientation('Euler',0,0,0,cs{c},ss);
%     end
%     
%     %Define Miller planes or directions
%     % if crys.Miller{c} is a cell array of either mtex miller objects, or numeric vectors
%     if isa(crys.Miller{c},'cell')                           
%         for i = 1:length(crys.Miller{c})                       %Loop over Miller indicee sets
%             if isa(crys.Miller{c}{i},'Miller')                 % mtex miller object - direct copy
%                 
%                 crys.oMil{c}{i} = crys.Miller{c}{i};
%                 crys.strMil{c}{i} = regexprep(num2str(round(...
%                     crys.oMil{c}{i}.(crys.type{c}))),'\s+',''); %Save Miller labels
%                 crys.Miller{c}{i}=crys.Miller{c}{i}.(crys.type{c});                      %overwrite crys.Miller to maintain compatibility 
%             elseif isa(crys.Miller{c}{i},'numeric') && length(crys.Miller{c}{i})==3   % 3-ix miller
%                 
%                 crys.oMil{c}{i} = Miller(crys.Miller{c}{i}(1),crys.Miller{c}{i}(2),...
%                     crys.Miller{c}{i}(3),cs{c},crys.type{c});  %Create Miller indexed direction
%                 crys.strMil{c}{i} = regexprep(num2str(round(...
%                     crys.oMil{c}{i}.(crys.type{c}))),'\s+',''); %Save Miller labels
%                 
%             elseif isa(crys.Miller{c}{i},'numeric') && length(crys.Miller{c}{i})==4   % 4-ix miller
%                 
%                 crys.oMil{c}{i} = Miller(crys.Miller{c}{i}(1),crys.Miller{c}{i}(2),...
%                     crys.Miller{c}{i}(3),crys.Miller{c}{i}(4),cs{c},crys.type{c});  %Create Miller indexed direction
%                 crys.strMil{c}{i} = regexprep(num2str(round(...
%                     crys.oMil{c}{i}.(crys.type{c}))),'\s+',''); %Save Miller labels
%             end
%             error(['Wrong input type for miller indices in crystal ' num2str(c) ', vector ' num2str(i)]);
%         end
%     % if crys.Miller{c} is a numeric vector of miller indices
%     elseif isa(crys.Miller{c},'numeric')                       
%         for i = 1:size(crys.Miller{c},1) %Loop over Miller indicee sets
%             if size(crys.Miller{c},2)==3 %3-index notation
%                 crys.oMil{c}{i} = Miller(crys.Miller{c}(i,1),crys.Miller{c}(i,2),...
%                     crys.Miller{c}(i,3),cs{c},crys.type{c});  %Create Miller indexed direction
%                 crys.strMil{c}{i} = regexprep(num2str(round(...
%                     crys.oMil{c}{i}.(crys.type{c}))),'\s+',''); %Save Miller labels
%             elseif  size(crys.Miller{c},2)==4 %4 index notation
%                 crys.oMil{c}{i} = Miller(crys.Miller{c}(i,1),crys.Miller{c}(i,2),...
%                     crys.Miller{c}(i,3),crys.Miller{c}(i,4),cs{c},crys.type{c});  %Create Miller indexed direction
%                 crys.strMil{c}{i} = regexprep(num2str(round(...
%                     crys.oMil{c}{i}.(crys.type{c}))),'\s+',''); %Save Miller labels
%             else
%                 error(['Wrong input type for miller indices in crystal ' num2str(c) ', vector ' num2str(i)]);
%             end
%         end
%     end
%     
% end

% *** Plot stereographic projection and inverse polefigure of equivalent crystal directions
% Reference frame is stage coordinate system C_s or equivalently microscope
% coordinate system C_m for 0ï¿½ stage tilt, no pre-tilt and no change in
% rotation with respect to the orientation measurement
if optim.plot
    for ii = 1:optim.order
        fprintf(['   -> Plotting stereographic projections and IPF of',...
            ' crystal %.0f ...\n'],ii);                                   %ScreenPrint
        stereoProj(crys{ii}.ori,crys{ii}.Miller,...
            sprintf('Original orientation %.0f',ii));                %Plot stereographic projection of equivalent crystal directions of crystal orientation 'o'
        figure;                                                            %Create figure
        plotIPDF(crys{ii}.ori,[xvector,yvector,zvector],'antipodal',...
            'MarkerFaceColor','k');                                   %Plot inverse pole-figure (IPF) [x:(100), y:(010) and z:(001)]
        set(gcf,'name',sprintf('Original orientation %.0f',ii));
    end
    tileFigs();                                                                %Sort figures
end
end

