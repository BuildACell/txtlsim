%FINDSPECIES  Find indices for species in a model or compartment
%
% findspecies(Mobj, namelist) returns a list of indices for the
% species given in namelist.  The names should be specified as a
% cell-array.

% RMM, 7 Sep 2012

function indexlist = findspecies(Mobj, namelist)

[nrows, len] = size(namelist);
indexlist = zeros(1, len);
for i=1:len
  for j =1:length(Mobj.Species)
    if strcmp(namelist(i), Mobj.Species(j).Name)
      indexlist(i) = j;
    end
  end
end
  