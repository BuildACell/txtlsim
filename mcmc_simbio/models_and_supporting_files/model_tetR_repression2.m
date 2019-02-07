function mobj = model_tetR_repression2
% repression with enzymatic one step protein production
% 
% ~~~ MODEL ~~~
% dT + pol <-> dT_pol -> dT + pol + mT
% mT + ribo <-> mT_ribo -> mT + ribo + pT
% mT -> null
%
% dG + pol <-> dG_pol -> dG + pol + mG
% mG + ribo <-> mG_ribo -> mG + ribo + pG
% mG -> null
% 
% 2 pT <-> pT2
% dG + pT2 <-> dG_pT2
% 
% Copyright (c) 2018, Vipul Singhal, Caltech
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

p = inputParser ;
addParameter(p, 'simtime', 2*3600);
parse(p);
p = p.Results;

%% setup model
mobj = sbiomodel('tetRepression2');

% Model Parameter values
cpol = 100; % nM
cribo = 50; %nM

rkfdG = 10; % nM-1s-1
rkrdG = 600; % s-1

rkcm = 0.001; %s-1

rkfdT = 10;
rkrdT = 600;

rkfpG = 10; % nM-1s-1
rkrpG = 300; % s-1

rkcp = 1/36;

rkfpT = 10;
rkrpT = 300;

rkfdimTet = 200; % nM-1s-1
rkrdimTet = 100; % s-1

rkfseqTet = 200; % nM-1s-1
rkrseqTet = 100; % s-1

rdel_m = log(2)/720; % 12 min half life of mrna

%% setup model reactions

% tetR TXTL
rxn = addreaction(mobj,'dT + pol <-> dT_pol');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfdT','krdT'};
addparameter(mobj, 'kfdT', rkfdT);
addparameter(mobj, 'krdT', rkrdT);

rxn = addreaction(mobj,'dT_pol -> dT + pol + mT');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kcm'};
addparameter(mobj, 'kcm', rkcm);

rxn = addreaction(mobj,'mT + ribo <-> mT_ribo');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfpT','krpT'};
addparameter(mobj, 'kfpT', rkfpT);
addparameter(mobj, 'krpT', rkrpT);

rxn = addreaction(mobj,'mT_ribo -> mT + ribo + pT');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kcp'};
addparameter(mobj, 'kcp', rkcp);

rxn = addreaction(mobj,'mT -> null');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'del_m'};
addparameter(mobj, 'del_m', rdel_m );

% GFP TXTL
rxn = addreaction(mobj,'dG + pol <-> dG_pol');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfdG','krdG'};
addparameter(mobj, 'kfdG', rkfdG);
addparameter(mobj, 'krdG', rkrdG);

rxn = addreaction(mobj,'dG_pol -> dG + pol + mG');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kcm'};

rxn = addreaction(mobj,'mG + ribo <-> mG_ribo');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfpG','krpG'};
addparameter(mobj, 'kfpG', rkfpG);
addparameter(mobj, 'krpG', rkrpG);

rxn = addreaction(mobj,'mG_ribo -> mG + ribo + pG');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kcp'};

rxn = addreaction(mobj,'mG -> null');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'del_m'};


rxn = addreaction(mobj,'2 pT <-> pT2');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfdimTet','krdimTet'};
addparameter(mobj, 'kfdimTet', rkfdimTet);
addparameter(mobj, 'krdimTet', rkrdimTet);

rxn = addreaction(mobj,'dG + pT2 <-> dG_pT2');
Kobj = addkineticlaw(rxn,'MassAction');
Kobj.ParameterVariableNames = {'kfseqTet','krseqTet'};
addparameter(mobj, 'kfseqTet', rkfseqTet);
addparameter(mobj, 'krseqTet', rkrseqTet);

% rxn = addreaction(mobj,'pT -> null');
% Kobj = addkineticlaw(rxn,'MassAction');
% Kobj.ParameterVariableNames = {'del_pT'};
% addparameter(mobj, 'del_pT', 0.1)

% rxn = addreaction(mobj,'pG -> null');
% Kobj = addkineticlaw(rxn,'MassAction');
% Kobj.ParameterVariableNames = {'del_pG'};
% addparameter(mobj, 'del_pG', 0.1)

% setup model species initial concentrations. 
% setup model species initial concentrations. 
specie = sbioselect(mobj, 'name', 'dT');
specie.InitialAmount = 1;

specie = sbioselect(mobj, 'name', 'dG');
specie.InitialAmount = 30;

specie = sbioselect(mobj, 'name', 'pT');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'pG');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'mT');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'mG');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'pol');
specie.InitialAmount = cpol;

specie = sbioselect(mobj, 'name', 'ribo');
specie.InitialAmount = cribo;

specie = sbioselect(mobj, 'name', 'dT_pol');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'dG_pol');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'mT_ribo');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'mG_ribo');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'pT2');
specie.InitialAmount = 0;

specie = sbioselect(mobj, 'name', 'dG_pT2');
specie.InitialAmount = 0;

%% Run the model

cs = getconfigset(mobj, 'active');
set(cs, 'StopTime', p.simtime);
        
sd = sbiosimulate(mobj);

end

