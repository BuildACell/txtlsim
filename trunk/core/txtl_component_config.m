classdef txtl_component_config
    %TXTL_COMPONTENT_CONFIG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        RNAPbound_Forward
        RNAPbound_Reverse
        Promoter_Length
        Thio_Length
        Junk_Length
        Gene_Length
        LVA_tag_Length
        Terminator_Length
        Dimmerization_Forward
        Dimmerization_Reverse
        Tetramerization_Forward
        Tetramerization_Reverse
        Protein_Inducer_Forward
        Protein_Inducer_Reverse
        Inducer_Degradation
        Protein_DNA_Forward
        Protein_DNA_Reverse
        Protein_Protein_Forward
        Protein_Protein_Reverse
        Protein_Maturation
    end
    
    methods
        function compConf = txtl_component_config(name)
            listOfProperties = properties(compConf);
            fileName = ['txtl_param_' name '.csv'];
            if exist(fileName,'file') == 2
                fid = fopen(fileName);
                % current parameter file format
                % Param_name, Param_type, Value, Comment
                M =  textscan(fid,'%q%q%q%q', 'delimiter',',','CollectOutput',true);
                M = M{1};
                fclose(fid);

                % set the object properties from the file.
                % Numeric/Expression type parameters are distinguished.
                for k = 1:size(listOfProperties)
                    index = find(cellfun(@(x) strcmp(x,listOfProperties(k)),M(:,1)) > 0);
                    if index > 0
                        if strcmp(M{index,2},'Numeric')
                            eval(sprintf('compConf.%s = %s;',M{index,1},M{index,3}));
                        else
                            % trying to evaluate expressions, if fails
                            % expression saved as string for later
                            % clarification
                            try
                               eval(sprintf('compConf.%s = %s;',M{index,1},M{index,3})); 
                            catch err
                               eval(sprintf('compConf.%s = ''%s'';',M{index,1},M{index,3}));
                            end
                        end
                    end
                end
            else
                error('the file: %s does not exist!',name);
            end
        end % end of constructor
        
    end
    
end

