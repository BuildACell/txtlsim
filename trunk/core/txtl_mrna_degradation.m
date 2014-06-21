function txtl_mrna_degradation(mode, tube, dna, rna, rbs_spec)
% Vipul Singhal Feb 13 2014
% We model RNA degradation as an enzymatic step reaction, as described in
% Vincent Noireaux's 2011 PRL paper 
% "Coarse-Grained Dynamics of Protein Synthesis in a Cell-Free System"
% 200 nM of radioactively labelled 960bases long RNA decause exponentially
% with a half life of 12 min. This means that the RNAse does not get
% saturated at this concentration of mRNA. 


% includeRNase = false;
% BooleanAtt1Present = false;
% BooleanAtt2Present = false;
% BooleanAnti1Present = false;
% BooleanAnti2Present = false;
% for i = 1: length(tube.userdata.DNAinfo)
%     UTR = tube.userdata.DNAinfo{i}{2};
%     [utrData, utrStr] = txtl_parsespec(UTR);
%     [BooleanAtt1Present,indexesOfAtt1Present] = checkForStringInACellList(utrData(1,:),'att1');
%     [BooleanAnti1Present,indexesOfAnti1Present] = checkForStringInACellList(utrData(1,:),'anti1');
%         [BooleanAtt2Present,indexesOfAtt2Present] = checkForStringInACellList(utrData(1,:),'att2');
%     [BooleanAnti2Present,indexesOfAnti2Present] = checkForStringInACellList(utrData(1,:),'anti2');
% if BooleanAtt1Present || BooleanAtt2Present || BooleanAnti1Present || BooleanAnti2Present
%     includeRNase = true;
% end
% end

includeRNase = true;

    complexF = tube.UserData.ReactionConfig.RNase_F;
    complexR = tube.UserData.ReactionConfig.RNase_R;
    degRate = tube.UserData.ReactionConfig.RNA_deg;
    if isempty(tube.UserData.ReactionConfig.RNase_F)
        complexF = degRate/10;
        complexR = degRate/40;
    end
    txtl_addreaction(tube,[rna.Name ' + RNase <-> ' rna.Name ':RNase'],...
        'MassAction',{'TXTL_RNAdeg_F',complexF;...
        'TXTL_RNAdeg_R',complexR});
    txtl_addreaction(tube,[rna.Name ':RNase -> RNase'],...
        'MassAction',{'TXTL_RNAdeg_F',degRate});
    
    
    if mode.utr_attenuator_flag
        % this needs to handle all the cases. So SIM, double antisense and
        % nominal. And later generalized string parsing. 
        [utrData, utrStr] = txtl_parsespec(rbs_spec);
        att = utrData{1,1};
        
        if mode.double_antisense
            %double antisense can only be anti1, no support for anti2 being
            %double yet. 
        txtl_addreaction(tube,['RNA ' att ' + RNase <-> RNA ' att ':RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',complexF;...
            'TXTL_RNAdeg_R',complexR});
        txtl_addreaction(tube,['RNA ' att ':RNase -> RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',degRate});
        
        txtl_addreaction(tube,['RNA ' att '-anti1 + RNase <-> RNA ' att '-anti1:RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',complexF;...
            'TXTL_RNAdeg_R',complexR});
        txtl_addreaction(tube,['RNA ' att '-anti1:RNase -> RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',degRate});
        txtl_addreaction(tube,['RNA ' att '-anti1-anti1 + RNase <-> RNA ' att '-anti1-anti1:RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',complexF;...
            'TXTL_RNAdeg_R',complexR});
        txtl_addreaction(tube,['RNA ' att '-anti1-anti1:RNase -> RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',degRate});
        elseif mode.sim_module_exception
        txtl_addreaction(tube,['RNA att1_SIM + RNase <-> RNA att1_SIM:RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',complexF;...
            'TXTL_RNAdeg_R',complexR});
        txtl_addreaction(tube,['RNA att1_SIM:RNase -> RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',degRate});    
        txtl_addreaction(tube,['RNA att1-att1 + RNase <-> RNA att1-att1:RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',complexF;...
            'TXTL_RNAdeg_R',complexR});
        txtl_addreaction(tube,['RNA att1-att1:RNase -> RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',degRate});     
        else
            
        txtl_addreaction(tube,['RNA ' att ' + RNase <-> RNA ' att ':RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',complexF;...
            'TXTL_RNAdeg_R',complexR});
        txtl_addreaction(tube,['RNA ' att ':RNase -> RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',degRate});
        end
    end
    
    
    if mode.utr_rbs_flag
        txtl_addreaction(tube,['AA:AGTP:Ribo:' rna.Name ' + RNase <-> AA:AGTP:Ribo:' rna.Name ':RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',complexF;...
            'TXTL_RNAdeg_R',complexR});
        txtl_addreaction(tube,['AA:AGTP:Ribo:' rna.Name ':RNase -> AGTP + AA + Ribo + RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',degRate});
        
        txtl_addreaction(tube,['Ribo:' rna.Name ' + RNase <-> Ribo:' rna.Name ':RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',complexF;...
            'TXTL_RNAdeg_R',complexR});
        txtl_addreaction(tube,['Ribo:' rna.Name ':RNase -> Ribo + RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',degRate});
    end

end

function [binVariable,indexes] = checkForStringInACellList(cellList,matchStr)
FlagVector = cellfun(@(x) strcmp(x,matchStr),cellList,'UniformOutput',false);
indexes = find(cell2mat(FlagVector) > 0);
if sum(cell2mat(FlagVector)) >= 1
    binVariable = true;
else
    binVariable = false;
end
end