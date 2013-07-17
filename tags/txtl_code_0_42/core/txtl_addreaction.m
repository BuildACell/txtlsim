function txtl_addreaction(tube,reactionEq,kineticLaw,parameters)


Robj1 = addreaction(tube, reactionEq);
Kobj1 = addkineticlaw(Robj1, kineticLaw);

for k=1:size(parameters,1)
    addparameter(Kobj1, parameters{k,1}, parameters{k,2});
end

set(Kobj1, 'ParameterVariableNames', parameters(:,1)');


end