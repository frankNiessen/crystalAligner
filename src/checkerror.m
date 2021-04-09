function crys = checkerror(crys)
assert(length(crys.alignAx) == length(crys.Miller), 'Number of axis-crystalDirection pairs of alignment objetives not equal');  %Check number of alignment axes and alignment objectives
if isa(crys.o,'cell')
    assert(length(crys.cs) == length(crys.o),'Number of CrystalSystem-Orientation pairs not equal');  %Check number of crystal systems and orientations
else
    assert(length(crys.cs) == size(crys.o,1),'Number of CrystalSystem-Orientation pairs not equal');  %Check number of crystal systems and orientations
end
%Check for different Miller notations
% safer to handle this natively in MTEX!!!
% edited defOri() to handle this
% %{
% for i = 1:length(crys.Miller)
%    assert(any(strcmpi(crys.type{i},{'hkl','hkil','uvw','uvtw'})),'Please choose ''uvw'', ''uvtw'', ''hkl'' or ''hkil'' as valid Miller type');
%    assert(size(crys.Miller{i},2) == 3 || size(crys.Miller{i},2) == 4,'Only 3 or 4 Miller indicees per set allowed');
%    if size(crys.Miller{i},2) == 4 && any(strcmpi(crys.type(i),{'hkl','hkil'})) %hkil notation
%       crys.Miller{i} = crys.Miller{i}(:,[1 2 4]);
%       crys.type{i} = 'hkl';
%    elseif size(crys.Miller{i},2) == 4 && any(strcmpi(crys.type(i),{'uvw','uvtw'})) %uvtw notation
%       crys.Miller{i} = [crys.Miller{i}(:,1)-crys.Miller{i}(:,3),crys.Miller{i}(:,2)-crys.Miller{i}(:,3),crys.Miller{i}(:,4)];
%       crys.type{i} = 'uvw';
%    end
% end
% %}
end