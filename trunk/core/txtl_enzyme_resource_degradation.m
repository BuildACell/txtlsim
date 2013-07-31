function txtl_enzyme_resource_degradation(modelObj)

% RNAP degradation in batch mode
if modelObj.UserData.Vesicule == 0
% RNAP degration as a first order reaction
RNAP_deg = 0.0011;
%
% % RNAP degradation
txtl_addreaction(modelObj,'RNAP -> null',...
    'MassAction',{'RNAPdeg_F',RNAP_deg});
txtl_addreaction(modelObj,'RNAP70 -> protein sigma70',...
    'MassAction',{'RNAPdeg_F',RNAP_deg});
txtl_addreaction(modelObj,'RNAP28 -> protein sigma28',...
    'MassAction',{'RNAPdeg_F',RNAP_deg});

end

%
% After 3hours because of the ATP regeneration stops the remaining ATP
% becomes unusable c.f. V Noireaux 2003.
atp_deg = 0.00003;

txtl_addreaction(modelObj,'ATP -> ATP_UNUSE',...
    'MassAction',{'ATPdeg_F',atp_deg});
% txtl_addreaction(modelObj,'AA:ATP:Ribo:RNA rbs--deGFP -> ATP_UNUSE + AA + Ribo:RNA rbs--deGFP',...
%     'MassAction',{'ATPdeg_F',atp_deg});
% txtl_addspecies(modelObj, 'ATP_REGEN_SUP',1);
%
% txtl_addreaction(modelObj,'ATP_REGEN_SUP -> null',...
%     'MassAction',{'ATP_F',0.00035});
%
% reactionObj = addreaction(modelObj, 'ATP_UNUSE -> ATP_REGEN_SUP + ATP');
% kineticlawObj = addkineticlaw(reactionObj, 'Henri-Michaelis-Menten');
%
%
% parameterObj1 = addparameter(kineticlawObj, 'Vm_d','Value',4);
% parameterObj2 = addparameter(kineticlawObj, 'Km_d','Value',1.25);
%
% set(kineticlawObj,'ParameterVariableNames', {'Vm_d' 'Km_d'});
% set(kineticlawObj,'SpeciesVariableNames', {'ATP_UNUSE'});

% txtl_addreaction(modelObj,'ATP_UNUSE:ATP_REGEN_SUP -> ATP_UNUSE',...
%     'MassAction',{'ATP_F',0.00035});
%
% txtl_addreaction(modelObj,'ATP_UNUSE + ATP_REGEN_SUP <-> ATP_UNUSE:ATP_REGEN_SUP',...
%     'MassAction',{'ATPregen_F',50; 'ATPregen_R',0.001});
%
% txtl_addreaction(modelObj,'ATP_UNUSE:ATP_REGEN_SUP -> ATP + ATP_REGEN_SUP',...
%     'MassAction',{'ATPregen_cat',30});

end