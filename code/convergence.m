close all
format shortE

func_num=0;
D=4;
Xmin=-2*pi;
Xmax=2*pi;
pop_size=100;
fhd = @fit_fun;
runs = 100;

iter_maxs=[40 80 120 160 200];

cur_date = datetime;
cur_date.Format = 'uuuu-MM-dd-HH-mm-ss';
cur_date = ['rast4d-' char(cur_date)];

dataq = parallel.pool.DataQueue;
% After receiving new data, update_progress() will be called
afterEach(dataq, @update_progress);

parfor k=1:length(iter_maxs)
    iter_max = iter_maxs(k);
    my_csv = zeros(runs,5);
    for i=1:runs
        tic
        [gbest,gbestval,FES, iter] = iLOA_func(fhd,D,pop_size,iter_max,Xmin,Xmax,[],cur_date,NaN,func_num);
        run_timed = toc;
        my_csv(i,4) = run_timed;
        my_csv(i,1) = gbestval;
        my_csv(i,2) = FES;
        my_csv(i,3) = iter;
        my_csv(i,5) = iter < iter_max;
        send(dataq,[iter_max i run_timed iter]);
    end
    csvwrite(['out/conv-test-' cur_date '-i_' num2str(iter_max) '.csv'], my_csv);
end

% Put this function at a proper location
function update_progress(i)
    if mod(i(2),10) == 0
        clc
    end
    fprintf('Finished: i_%d - %d at %g iters %d\n', i(1), i(2), i(3), i(4));
end

% function fitness = fit_fun(pos, n)
%     fitness = (pos(1)*pi/100)^2;
% end

% RASTRIGIN (0, 0) minima n(5) [0,50]
function fitness = fit_fun(pos, n)
    dimensions = length(pos);
    fitness = 10*dimensions;
    for i=1:dimensions
        xi = pos(i);
        fitness = fitness + xi ^ 2 - 10 * cos(2*pi*xi);
    end
end

% ROSENBROCK2d (0, 0) minima n(5) [0,10000]
% function fitness = fit_fun(pos, n)
%     fitness = (1-pos(1))^2 + 100*(pos(2) - pos(1)^2)^2;
% end

% function fitness = fit_fun(pos, n)
%     dimensions = length(pos);
%     fitness = 0;
%     for i=1:dimensions
%         fitness = fitness + pos(i)^2;
%     end
% end

% Griewank 1d*10 n(100) [0,5]
% function fitness = fit_fun(pos, n)
%     fitness = 1 + (1/4000)*pos(1)^2-cos(pos(1)*10);
% end

% Griewank 2d
% function fitness = fit_fun(pos, n)
%     fitness = 1 + (1/4000)*pos(1)^2 + (1/4000)*pos(2)^2-cos(pos(1)) * cos(sqrt(2)*pos(2)/2);
% end

% Griewank 3d
% function fitness = fit_fun(pos, n)
%     fitness = 1 + (1/4000)*pos(1)^2 + (1/4000)*pos(2)^2 + (1/4000)*pos(3)^2 - cos(pos(1)) * cos(sqrt(2)*pos(2)/2) * cos(3^(1/3) * pos(3) / 3);
% end