%FINDSPECIES  Find indices for species in a model or compartment
%
% findspecies(Mobj, namelist) returns a list of indices for the
% species given in namelist.  The names should be specified as a
% cell-array.

% RMM, 7 Sep 2012

function indexlist = findspecies(Mobj, namelist)

if isstr(namelist)
  % Initialize the return value to zero, in case we don't find anything
  indexlist = 0;

  
  %! TODO zoltuz 1/10/12 cellfun for speed up: find(cellfun(@(x)
  %strcmp('AA',x),listOfSpecies) == 1)
  
  % Search through the list for our species name
  for j =1:length(Mobj.Species)
    if strcmp(namelist, Mobj.Species(j).Name)
      indexlist = j;
      return
    end
  end

else
    
  % !TODO zoltuz 1/10/12 intersect for speed up  
  % Return an array of data of same size as input
  [nrows, len] = size(namelist);
  indexlist = zeros(1, len);
  for i=1:len
    for j =1:length(Mobj.Species)
      if strcmp(namelist(i), Mobj.Species(j).Name)
	indexlist(i) = j;
      end
    end
  end
end

% Local variables:
% mode: matlab
% End:
