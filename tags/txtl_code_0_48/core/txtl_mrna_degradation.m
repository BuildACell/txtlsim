function txtl_mrna_degradation(mode, tube, dna, rna, rbs_spec)


% since RNase is only used in the JLucks RNA circuits, figure out if an att
% or anti is present in ANY of the DNAs. if so, then set up degradations using RNase. 
includeRNase = false;
BooleanAtt1Present = false;
BooleanAtt2Present = false;
BooleanAnti1Present = false;
BooleanAnti2Present = false;
for i = 1: length(tube.userdata.DNAinfo)
    UTR = tube.userdata.DNAinfo{i}{2};
    [utrData, utrStr] = txtl_parsespec(UTR);
    [BooleanAtt1Present,indexesOfAtt1Present] = checkForStringInACellList(utrData(1,:),'att1');
    [BooleanAnti1Present,indexesOfAnti1Present] = checkForStringInACellList(utrData(1,:),'anti1');
        [BooleanAtt2Present,indexesOfAtt2Present] = checkForStringInACellList(utrData(1,:),'att2');
    [BooleanAnti2Present,indexesOfAnti2Present] = checkForStringInACellList(utrData(1,:),'anti2');
if BooleanAtt1Present || BooleanAtt2Present || BooleanAnti1Present || BooleanAnti2Present
    includeRNase = true;
end
end



if includeRNase
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
        [utrData, utrStr] = txtl_parsespec(rbs_spec);
        att = utrData{1,1};
        txtl_addreaction(tube,['RNA ' att ' + RNase <-> RNA ' att ':RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',complexF;...
            'TXTL_RNAdeg_R',complexR});
        txtl_addreaction(tube,['RNA ' att ':RNase -> RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',degRate});
        if mode.sim_module_exception
            txtl_addreaction(tube,'RNA att1-att1 + RNase <-> RNA att1-att1:RNase',...
                'MassAction',{'TXTL_RNAdeg_F',complexF;...
                'TXTL_RNAdeg_R',complexR});
            txtl_addreaction(tube,'RNA att1-att1:RNase -> RNase',...
                'MassAction',{'TXTL_RNAdeg_F',degRate});
        elseif mode.double_antisense
            txtl_addreaction(tube,'RNA att2-anti1 + RNase <-> RNA att2-anti1:RNase',...
                'MassAction',{'TXTL_RNAdeg_F',complexF;...
                'TXTL_RNAdeg_R',complexR});
            txtl_addreaction(tube,'RNA att2-anti1:RNase -> RNase',...
                'MassAction',{'TXTL_RNAdeg_F',degRate});
        end
    end
    
    
    if mode.utr_rbs_flag
        txtl_addreaction(tube,['AA:ATP:Ribo:' rna.Name ' + RNase <-> AA:ATP:Ribo:' rna.Name ':RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',complexF;...
            'TXTL_RNAdeg_R',complexR});
        txtl_addreaction(tube,['AA:ATP:Ribo:' rna.Name ':RNase -> ATP + AA + Ribo + RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',degRate});
        
        txtl_addreaction(tube,['Ribo:' rna.Name ' + RNase <-> Ribo:' rna.Name ':RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',complexF;...
            'TXTL_RNAdeg_R',complexR});
        txtl_addreaction(tube,['Ribo:' rna.Name ':RNase -> Ribo + RNase'],...
            'MassAction',{'TXTL_RNAdeg_F',degRate});
    end
else
        txtl_addreaction(tube,[rna.Name ' -> null'],...
        'MassAction',{'TXTL_RNAdeg_F',tube.UserData.ReactionConfig.RNA_deg});
        txtl_addreaction(tube,['AA:ATP:Ribo:' rna.Name ' -> AA + ATP + Ribo'],...
            'MassAction',{'TXTL_AA_ATP_RNAdeg_F',tube.UserData.ReactionConfig.RNA_deg});
        txtl_addreaction(tube,['Ribo:' rna.Name ' -> Ribo'],...
            'MassAction',{'TXTL_Ribo_RNAdeg_F',tube.UserData.ReactionConfig.RNA_deg});
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