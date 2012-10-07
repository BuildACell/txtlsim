function txtl_plot(t_ode,x_ode,modelObj,dataGroups)
% initial version for txtl_plot, the RNAs and proteins are automatically exploited from the provided DNA sequence. 
% t_ode: nx1 time vector, no time scaling is applied inside!
% x_ode: nxm species vector
% modelObj: simBiology object of the current model
% dataGroups: special data structure for the plots (regular expressions can added)
%
% example of the required data structure
%  * First column is the name of desired plot (available plots are listed
%    below)
%
%  * Second column contains name-list of the simulation data, which will be plotted. 
%    (RNA and protein names are exploited and plotted automatically from DNA sequences)  
%  +-----------------------------------------
%  |Special strings:
%  |
%  | * it handles matlab compatible regular expressions 
%  |   (e.g. plotting all of the protein in the system: dataGroups{2,2} = {'#(protein \w*)'};)
%  | * txtl_plot also can calculate the total concentration of selected
%  | proteins and its variants with sprint "[protein name]_tot", where protein
%  | is a valid Species name in the modelObj. (e.g. dataGroups{2,2} = {'[protein lacI]_tot'} )
%  |
%  +-----------------------------------------
%
%  * Third column is optional and it is designated for user defined line coloring
%
% Currently 3 types of plots are supported:
% - DNA and mRNA plot (case sensitive!)
% - Gene Expression plot
% - Resource usage
%
% first the DNA sequences should be provided for automatic name extraction
% %DNA and mRNA plot
%  dataGroups{1,1} = 'DNA and mRNA';
%  dataGroups{1,2} = {'DNA p70=rbs=lacI','DNA placi=rbs=deGFP'}%,'RNA rbs=lacI','RNA rbs=deGFP'}
%  dataGroups{1,3} = {'b-','r-','b--','r--'}
%
%
%
% %Gene Expression Plot
%  dataGroups{2,1} = 'Gene Expression';
%  dataGroups{2,2} = {'protein deGFP*','protein gamS','protein lacIdimer', 'protein lacItetramer'};
%  dataGroups{2,3} = {'b-','g--','g-','r-','b--','b-.'}
%
% 
% %Resource Plot
%  dataGroups{3,1} = 'Resource usage';
%

%%


numOfGroups = size(dataGroups,1);
listOfProteins = {};
listOfRNAs = {};
listOfDNAs = {};

figure(1); clf(); 

for k = 1:numOfGroups

   %%%%%%%
   % DNA and mRNA plot
   %
   %%%%%%%
   if(strcmp(dataGroups{k,1},'DNA and mRNA'))

    %! TODO further refinement of str spliting
    if ~isempty(dataGroups{k,2})
        regexp_ind = strmatch('#', dataGroups{k,2});
        if ~isempty(regexp_ind)
            autoSpecies = extractRegexpAndExecute(regexp_ind,modelObj,dataGroups{k,2});
            % remove regexp str form dataGroups array
            dataGroups{k,2}(regexp_ind) = [];
            listOfDNAs =  horzcat(listOfDNAs,autoSpecies);
            % delete multiple occurrences
            listOfDNAs =  unique(listOfDNAs);
        end
        listOfDNAs = horzcat(listOfDNAs,dataGroups{k,2});
        
        r = regexp(listOfDNAs,'=','split');
        %! TODO skip already added proteins and mRNAs to avoid duplicates
        % get each RNAs
        RNAs = cellfun(@(x) strcat('RNA',{' '},x(2),'=',x(3)),r);
        % get each proteins
        proteins = cellfun(@(x) strcat('protein',{' '}, x(3)), r); % adding protein string for each element
        listOfProteins = horzcat(listOfProteins,proteins);
        listOfRNAs = horzcat(listOfRNAs,RNAs);

        % collect the data   
        listOfDNAsRNAs = horzcat(listOfDNAs,listOfRNAs);
        dataX = getDataForSpecies(modelObj,x_ode,listOfDNAsRNAs);
    else
        warning('No DNA strings were provided!');
    end
    % plot the data 
     subplot(223)
    if (~isempty(dataGroups{k,3}))
      [ColorMtx,LineStyle] = getColorAndLineOrderByUserData(dataGroups{k,3});
      for l=1:size(dataX,2)
          line('XData',t_ode/60,'YData',dataX(:,l),'Color',ColorMtx(l,:),'LineStyle',LineStyle{l});
      end
    else
      plotID_dna = plot(t_ode/60,dataX);
    end
    
   
    
    lgh = legend(listOfDNAsRNAs, 'Location', 'Best');
    legend(lgh, 'boxoff');
    ylabel('Species amounts [nM]');
    xlabel('Time [min]');
    title(dataGroups{k,1});
   %%%%%%%%%%
   % Gene Expression plot
   %
   %%%%%%%%%%
   elseif(strcmp(dataGroups{k,1},'Gene Expression'))
       
    % if total amount of selected protein is calculated   
    totalAmount = {};  
    % add extra user defined proteins
    if ~isempty(dataGroups{k,2}) 
        regexp_ind = strmatch('#', dataGroups{k,2});
        matchStrings = regexp(dataGroups{k,2},'^\[(protein \w*)\]_tot','tokens')
        if ~isempty(regexp_ind)
           % if regular expression is present, than it is executed and the result goes to the listOfProteins 
           auto_species =  extractRegexpAndExecute(regexp_ind,modelObj,dataGroups{k,2});
           % remove regexp str form dataGroups array
           dataGroups{k,2}(regexp_ind) = [];
           
           listOfProteins =  horzcat(listOfProteins,dataGroups{k,2},auto_species);
           % delete multiple occurrences
           listOfProteins =  unique(listOfProteins);
           
        % calculating the total concentration of selected proteins   
        elseif ~isempty(~cellfun(@(x) isempty(x),matchStrings))
            for z = 1:size(matchStrings,2)
              u = 1;  
              if ~isempty(matchStrings{z})
                dataGroups{k,2}(z) = [];  
                totalAmount{u,1} = matchStrings{z}{1}{1};
                totalAmount{u,2} = [];
               %! TODO is this the best way to do that? 
               for l = 1:size(modelObj.Species,1)
                   regString = [matchStrings{z}{1}{1} '.*'];
                   matchProtein = regexp(modelObj.Species(l).Name,regString,'match');
                   if ~isempty(matchProtein)
                       if isempty(totalAmount{u,2})
                          totalAmount{u,2} =  x_ode(:,l);
                       else
                          totalAmount{u,2} = totalAmount{u,2} + x_ode(:,l);
                       end
                   end
               end
               u = u +1;
              end
            end
        % no special string in the list, so we add it to the rest.    
        else
           listOfProteins =  horzcat(listOfProteins,dataGroups{k,2});
        end
       
        
    end
    dataX = getDataForSpecies(modelObj,x_ode,listOfProteins);
    % adding total protein concentraion into the common data matrix and a
    % label matrix as well (This could be done before, because the
    % listOfProteins was used for aquiring Species data by name)
    if size(totalAmount,1) > 0
        for k = 1:size(totalAmount,1)
        listOfProteins(end+1) = {sprintf('[%s]_{tot}',totalAmount{k,1})};     
        dataX(:,end+1) = totalAmount{k,2};
        end
    end 
   
 
   subplot(2,2,1:2);
    
        if (~isempty(dataGroups{k,3}))
            [ColorMtx,LineStyle] = getColorAndLineOrderByUserData(dataGroups{k,3});
            for l=1:size(dataX,2)
              line('XData',t_ode/60,'YData',dataX(:,l),'Color',ColorMtx(l,:),'LineStyle',LineStyle{l});
            end
        else
        plotID = plot(t_ode/60,dataX);
        end
    
   
    lgh = legend(listOfProteins, 'Location', 'NorthEast');
    legend(lgh, 'boxoff');
    ylabel('Species amounts [nM]');
    xlabel('Time [min]');
    title(dataGroups{k,1});
   
   %%%%%%%%%%
   % Resource usage plot
   %
   %%%%%%%%%% 
   elseif(strcmp(dataGroups{k,1},'Resource usage'))
       
    listOfResources = {'NTP','AA','RNAP','Ribo'};
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
      {'NTP [mM]', 'AA [mM]', 'RNAP [nM]', 'Ribo [nM]'}, ...
      'Location', 'Best');
    legend(lgh, 'boxoff');
    ylabel('Species amounts [normalized]');
    xlabel('Time [min]');

   else
       disp('no option were provided!');
   
   end % end of if dataGroups
    
end % end of for
    
end

function autoSpecies = extractRegexpAndExecute(regexp_ind,modelObj,dataSource)
autoSpecies = {};
    for z=1:size(regexp_ind,1)
        r_str = strrep(dataSource(z),'#','');
        for l=1:size(modelObj.Species,1)
          specie_match = regexp(modelObj.Species(l).Name,r_str,'tokens');
           if(~cellfun(@(x) isempty(x),specie_match))
           autoSpecies(end+1) = specie_match{1}{1};
           end
        end
    end       
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


