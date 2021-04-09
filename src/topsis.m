function ind = topsis(epsIn,w,normMode)
%function ind = topsis(epsIn,w)
%TOPSIS - Technique for Order of Preference by Similarity to Ideal Solution (TOPSIS)
%Multi-objective decision-making method to rank Pareto Solutions
%INPUT
%epsIn: Fitness function values [NrParetoSolutions x NrObjectives]
%w:     Objective weights [1 x NrObjectives]
%ind:   Indices of Pareto solutions sorted from highest to lowest score
if strcmpi(normMode,'Vector')
   epsNorm = epsIn./sqrt(sum(epsIn.^2));                                   %Normalized fitness function values by column
elseif strcmpi(normMode,'Matrix')
   epsNorm = epsIn./sqrt(sum(epsIn(:).^2));                                %Normalized fitness function values by matrix
else
   error(sprintf('Invalid normalization mode ''%s''',normMode));           %Error
end

M = epsNorm.*repmat(w,size(epsNorm,1),1);                                  %Compute decision matrix
Mworst = max(M,[],1);                                                      %Find maxima
Mbest = min(M,[],1);                                                       %Find minima
Lworst = sqrt(sum((M-repmat(Mworst,size(M,1),1))'.^2));                    %Compute L-square-distance between M and Mworst
Lbest = sqrt(sum((M-repmat(Mbest,size(M,1),1))'.^2));                      %Compute L-square-distance between M and Mbest
score = Lworst./(Lworst+Lbest);                                            %Compute score
[~,ind] = sort(score,'descend');                                           %Rank score
end
