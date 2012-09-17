function plotID = txtl_genexpression_plot(simObj,modelObj,listOfProtein)

if(isempty(listOfProtein))
   error('you must provide valid list of protein'); 
end

numOfProtein = size(listOfProtein,2);

indexNum = zeros(1,numOfProtein);
dataX = zeros(size(simObj.Time,1),numOfProtein);
for k = 1:numOfProtein
    indexNum(k) = findspecies(modelObj, ['protein ' listOfProtein{k}])
    if(indexNum(k) == 0)
        error('not valid protein Name: protein %s',listOfProtein{k});
    end
    dataX(:,k) = simObj.Data(:,indexNum(k));
end

%! TODO add color and line style:
% http://blogs.mathworks.com/loren/2007/12/19/plotting-with-style/

%! TODO just return with the ID and plot it latter on.
plotID = plot(simObj.Time/60,dataX)


title('Gene Expression');
lgh = legend(listOfProtein, 'Location', 'Best');
legend(lgh, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');



end