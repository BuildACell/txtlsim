% single parameter estimation for LacOperon


% First we want to estimate the production rate of allolactose

% resample time and data
resample_interval = 10;

[t_alloLac,x_alloLac] = selectbyname(simData, 'alloLactose');

t_sample = t_alloLac(1:resample_interval:end); 
x_sample = x_alloLac(1:resample_interval:end);

% add some noise
x_exp = x_sample + 0.01*randn(size(x_sample));

% compare orig and resampled
figure(1)
plot(t_alloLac/60,x_alloLac)
hold on
plot(t_sample/60,x_sample,'ro')
plot(t_sample/60,x_exp,'c*')



paramObj = sbioselect(well_a1, 'Name', 'Vl_alloLac')


 opt = optimset('PlotFcns',@optimplotfval,'MaxIter',15);
 [estValues1, result1] = sbioparamestim(well_a1, t_sample, x_exp, ...
     well_a1.Species(19),paramObj, {}, {'fminsearch',opt});
 
 
 
 % multiple parameter estimation