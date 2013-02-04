classdef txtl_reaction_config
    %TXTL_REACTION_CONFIG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        NTPmodel;
        AAmodel;
        Transcription_Rate;
        Translation_Rate;
        DNA_RecBCD_Forward;
        DNA_RecBCD_Reverse;
        DNA_RecBCD_complex_deg;
        Protein_ClpXP_Forward;
        Protein_ClpXP_Reverse;
        Protein_ClpXP_complex_deg;
        RNAP_S70_F;
        RNAP_S70_R;
        AA_Forward;
        AA_Reverse;
        Ribosome_Binding_F;
        Ribosome_Binding_R;
        RNA_deg;
        NTP_Forward;
        NTP_Reverse;
    end
    
    methods
        function rConf = txtl_reaction_config(name)
            listOfProperties = properties(rConf);
            fileName = [name '_config.csv'];
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
                            eval(sprintf('rConf.%s = %s;',M{index,1},M{index,3}));
                        else
                            % trying to evaluate expressions, if fails
                            % expression saved as string for later
                            % clarification
                            try
                               eval(sprintf('rConf.%s = %s;',M{index,1},M{index,3})); 
                            catch err
                               eval(sprintf('rConf.%s = ''%s'';',M{index,1},M{index,3}));
                            end
                        end
                    end
                end
            else
                error('the file: %s does not exist!',name);
            end
        end
        
        
      
    end
    
end

