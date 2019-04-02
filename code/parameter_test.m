close all

func_num=0;
D=2;
Xmin=-100;
Xmax=100;
pop_size=50;
iter_max=30;
fhd = @fit_fun;
runs = 5;

rand_lio = rand(D, pop_size).*(Xmax-Xmin)+Xmin;

cur_date = datetime;
cur_date.Format = 'uuuu-MM-dd-HH-mm-ss';
cur_date = char(cur_date);

dataq = parallel.pool.DataQueue;
% After receiving new data, update_progress() will be called
afterEach(dataq, @update_progress);

percentNomad = [0.2 0.4 0.6 0.8];
percentRoam = [0.2 0.4 0.6 0.8];
percentSex = [0.2 0.4 0.6 0.8];
rateMating = [0.3 0.6];
probabilityMutation = [0.2 0.4 0.6 0.8];
rateImmigration = [0.2 0.4 0.6 0.8];
percentGroupInfluence = [0.2 0.4 0.6 0.8];
pressureRankedSelect = [2 3];
pressureNearBest = [2 3];

rundata = [];

n_completed = 0;

parfor i=1:length(percentNomad)
    count = 0;
    cur_set = iLOA_setting;
    cur_set.numberOfPrides = 4;
    cur_set.annealing = true;
    cur_set.percentNomad = percentNomad(i);
    for j=1:length(percentRoam)
        cur_set.percentRoam = percentRoam(j);
        for k=1:length(percentSex)
            cur_set.percentSex = percentSex(k);
            for l=1:length(rateMating)
                cur_set.rateMating = rateMating(l);
                for r=1:length(probabilityMutation)
                    cur_set.probabilityMutation = probabilityMutation(r);
                    for m=1:length(rateImmigration)
                        cur_set.rateImmigration = rateImmigration(m);
                        for n=1:length(percentGroupInfluence)
                            cur_set.percentGroupInfluence = percentGroupInfluence(n);
                            for o=1:length(pressureRankedSelect)
                                cur_set.pressureRankedSelect = pressureRankedSelect(o);
                                for p=1:length(pressureNearBest)
                                    cur_set.pressureNearBest = pressureNearBest(p);
                                    runline = zeros(1,runs);
                                    for q=1:runs
                                        [gbest,gbestval,FES] = iLOA_func(fhd,D,pop_size,iter_max,Xmin,Xmax,rand_lio,cur_date,cur_set,func_num);
                                        runline(q) = gbestval;
                                    end
                                    runline = [cur_set.get_all() runline];
                                    count = count + 1;
                                    dlmwrite(['out/paramtest-' cur_date '.csv'],runline,'-append','delimiter',',','roffset',0,'precision',10)
                                    send(dataq,[i count]);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


% Put this function at a proper location
function update_progress(i)
    if mod(i(2),10) == 0
        clc
    end
    fprintf('Finished: %d - %d\n', i(1), i(2));
end

% function fitness = fit_fun(pos, n)
%     fitness = (pos(1)*pi/100)^2;
% end

% RASTRIGIN (0, 0) minima n(5) [0,50]
% function fitness = fit_fun(pos, n)
%     dimensions = length(pos);
%     fitness = 10*dimensions;
%     for i=1:dimensions
%         xi = pos(i);
%         fitness = fitness + xi ^ 2 - 10 * cos(2*pi*xi);
%     end
% end

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

% Griewank 1d n(100) [0,5]
% function fitness = fit_fun(pos, n)
%     fitness = 1 + (1/4000)*pos(1)^2-cos(pos(1));
% end

% Griewank 2d
function fitness = fit_fun(pos, n)
    fitness = 1 + (1/4000)*pos(1)^2 + (1/4000)*pos(2)^2-cos(pos(1)) * cos(sqrt(2)*pos(2)/2);
end

