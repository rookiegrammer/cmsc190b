clear
clc
% mex cec14_func.cpp -DWINDOWS
t = datetime('now','TimeZone','local','Format','d-MMM-y-HH_mm_ss');
tstr = sprintf('%s', t);
title = 'LOA';

func_num=1;
D=2;
Xmin=-100;
Xmax=100;
pop_size=50;
iter_max=50;
runs=1;
fhd=str2func('cec14_func');
for i=16:30 % functions
    func_num=i;
    fprintf('Function %d\nRunning %d...\n', i, runs);
    timer=tic;
    fprintf('\n');
    for j=1:runs
        [gbest,gbestval,FES]= LOA_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
        xbest(j, :)=gbest;
        fbest(i,j)=gbestval;
        fprintf('%d..', j);
    end
    f_min = min(fbest(i,:));
    file_id = fopen([title '-result-' tstr '-fun-' num2str(i) '.txt'], 'wt');
    fprintf(file_id, '%s', t);
    fprintf(file_id, '\n - Ran %f sec, min %g\n\n', toc(timer), f_min);
    fprintf(file_id, '\nBest Xs\n');
    for ii = 1:size(xbest, 1)
        fprintf(file_id, 'Run %02d: ',ii);
        fprintf(file_id, '%g\t', xbest(ii,:));
        fprintf(file_id, '\n');
    end
    fprintf(file_id, '\nBest Xs Fitness\n');
    fprintf(file_id, '%g\n', fbest(i, :));
    fclose(file_id);
end
fprintf('Done!');

% for i=1:30
% eval(['load input_data/shift_data_' num2str(i) '.txt']);
% eval(['O=shift_data_' num2str(i) '(1:10);']);
% f(i)=cec14_func(O',i);
% end
% f