function indexList = findStringInAList (list,string)

 if ~isempty(list) && ~isempty(string)
    if ischar(list)
        indexList = strcmp(list,string);
    else
        indexes = cellfun(@(x) strcmp(x,string),list);
        ind = find(indexes > 0);
        if ~isempty(ind)
            indexList = ind;
        else
            indexList = 0;
        end
    end
 else
     indexList = 0;
 end
end