function crys = checkerror(crys)
%function crys = checkerror(crys)
%Short error check

if length(crys)  == 1 || length(crys)  == 2
    assert(isfield(crys{1},'alignAx') &&  isfield(crys{1},'Miller') && length(crys{1}.alignAx)==1 && ~isempty(crys{1}.Miller),...
           'The first alignment objective needs a single microscope axis and at least one crystal directions');  
else
    error("You may only specify one or two alignment objectives!");
end
if length(crys)  == 2
    assert(isfield(crys{2},'alignAx') &&  isfield(crys{2},'Miller') && length(crys{2}.alignAx)==1 && length(crys{2}.Miller)==1,...
           'The second alignment objective needs a single microscope axis and one crystal directions');  
end


end
