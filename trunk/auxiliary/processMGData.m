function [whole_exp_file] = processMGData(whole_exp_file)

% search for mg channel
channel_index = cellfun(@(x) strcmpi(x,'MG'),whole_exp_file.channels(:,1)) > 0;

t_vec = whole_exp_file.t_vec/60;

disp(sprintf('calculating the MG curve fit for individual wells in %s...',whole_exp_file.FileName));
tic
f2 = fittype('gauss3');
for k = 1:size(whole_exp_file.Data,2);
    [fitinfo2(k).ff gof2(k).gof] = fit(t_vec,whole_exp_file.noBg(:,k,channel_index),f2);    
    whole_exp_file.MgCurve(:,k) = fitinfo2(k).ff(whole_exp_file.t_vec/60);
end
toc
whole_exp_file.MgCurveFit = fitinfo2;
whole_exp_file.MgCurveFitGOF = gof2;




end