function txtl_plot(t_ode,x_ode,modelObj,dataGroups)
% initial version for txtl_plot
% t_ode: nx1 time vector, no time scaling is applied inside!
% x_ode: nxm species vector
% modelObj:
% name: list of DNA sequnces to display

numOfGroups = size(dataGroups,1);
listOfProteins = {};
listOfRNAs = {};
listOfDNAs = {};

figure(1); clf(); 

for k = 1:numOfGroups

   if(strcmp(dataGroups{k,1},'DNA and mRNA'))

    %! TODO further refinement of str spliting
    r = regexp(dataGroups{k,2},'=','split');
    %! TODO skip already added proteins and mRNAs to avoid duplicates
    % get each RNAs
    RNAs = cellfun(@(x) strcat('RNA',{' '},x(2),'=',x(3)),r);
    % get each proteins
    proteins = cellfun(@(x) strcat('protein',{' '}, x(3)), r); % adding protein string for each element
    listOfProteins = horzcat(listOfProteins,proteins);
    listOfRNAs = horzcat(listOfRNAs,RNAs);
    
    % collect the data   
    listOfDNAs = dataGroups{k,2};
    listOfDNAsRNAs = horzcat(listOfDNAs,listOfRNAs)
    dataX = getDataForSpecies(modelObj,x_ode,listOfDNAsRNAs);

    % plot the data 
     subplot(223)
    if (~isempty(dataGroups{k,3}))
      [ColorMtx,LineStyle] = getColorAndLineOrderByUserData(dataGroups{k,3})
      
      for l=1:size(dataX,2)
            h(l) = line('XData',t_ode/60,'YData',dataX(:,l),'Color',ColorMtx(l,:),'LineStyle',LineStyle{l})
      end
      set(gca,'Children',[h])    
    else
      plotID_dna = plot(t_ode/60,dataX);
    end
    
   
    
    lgh = legend(listOfDNAsRNAs, 'Location', 'Best');
    legend(lgh, 'boxoff');
    ylabel('Species amounts [nM]');
    xlabel('Time [min]');
    title(dataGroups{k,1});
    
   elseif(strcmp(dataGroups{k,1},'Gene Expression'))
      
   % add extra user defined proteins
    listOfProteins =  horzcat(listOfProteins,dataGroups{k,2});
    dataX = getDataForSpecies(modelObj,x_ode,listOfProteins);
   
 
    subplot(221);
    
        if (~isempty(dataGroups{k,3}))
        %gObj = axes();   
        [ColorMtx,LineStyle] = getColorAndLineOrderByUserData(dataGroups{k,3})
        for l=1:size(dataX,2)
            h(l) = line('XData',t_ode/60,'YData',dataX(:,l),'Color',ColorMtx(l,:),'LineStyle',LineStyle{l})
        end
        set(gca,'Children',[ fliplr(h)])    
        else
        plotID = plot(t_ode/60,dataX);
        end
    
   
    lgh = legend(listOfProteins, 'Location', 'NorthWest');
    legend(lgh, 'boxoff');
    ylabel('Species amounts [nM]');
    xlabel('Time [min]');
    title(dataGroups{k,1});
   
    
   elseif(strcmp(dataGroups{k,1},'Resource usage'))
       
    listOfResources = {'NTP','AA','RNAP70','Ribo'};
    dataX = getDataForSpecies(modelObj,x_ode,listOfResources);
    
    subplot(224)
    mMperunit = 100 / 1000;			% convert from NTP, AA units to mM
    plot(...
      t_ode/60, dataX(:, 1)/dataX(1, 1), 'b-', ...
      t_ode/60, dataX(:, 2)/dataX(1, 2), 'r-', ...
      t_ode/60, dataX(:, 3)/dataX(1, 3), 'b--', ...
      t_ode/60, dataX(:, 4)/dataX(1, 4), 'r--');

    title('Resource usage');
    lgh = legend(...
      {'NTP [mM]', 'AA [mM]', 'RNAP70 [nM]', 'Ribo [nM]'}, ...
      'Location', 'Best');
    legend(lgh, 'boxoff');
    ylabel('Species amounts [normalized]');
    xlabel('Time [min]');

   else
       
   
   end % end of if dataGroups
    
end % end of for
    

end

function dataX = getDataForSpecies(modelObj,x_ode,listOfSpecies)
% collect data for the listed species from the simulation result array (x_ode)

  indexNum = findspecies(modelObj, listOfSpecies);
  notASpecie = find(indexNum == 0);
  if (~isempty(notASpecie))
     for k=1:size(notASpecie,2)
       error('not valid specie name: %s!',listOfSpecies{notASpecie(k)});
     end
  end
  dataX = x_ode(:,indexNum);
  
end

function rgbVector = convertLabelToRGBValue(label)
rgbVector = rem(floor((strfind('kbgcrmyw', label) - 1) * [0.25 0.5 1]), 2);

end

%! TODO handle only color expressions
function [ColorMtx,LineStyle] = getColorAndLineOrderByUserData(listOfItems)
        styleOptions = regexp(listOfItems,'(.)(.*)','tokens');
        ColorMtx = [];
        LineStyle = {};
  
        for l=1:size(styleOptions,2)
         ColorMtx(l,:) = convertLabelToRGBValue(styleOptions{1,l}{1,1}{1});
         LineStyle{l} = styleOptions{1,l}{1,1}{2}; 
        end
end

