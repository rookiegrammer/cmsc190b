close all

func_num=0;
D=2;
Xmin=-5;
Xmax=5;
pop_size=50;
iter_max=50;
fhd = @fit_fun;

rand_lio = rand(D, pop_size).*(Xmax-Xmin)+Xmin;

cur_date = datetime;
cur_date.Format = 'uuuu-MM-dd-HH-mm-ss';
cur_date = char(cur_date);

[gbest,gbestval,FES] = LOA_func(fhd,D,pop_size,iter_max,Xmin,Xmax,rand_lio,cur_date,func_num)
[gbest2,gbestval2,FES2] = iLOA_func(fhd,D,pop_size,iter_max,Xmin,Xmax,rand_lio,cur_date,NaN,func_num)

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
% function fitness = fit_fun(pos, n)
%     fitness = 1 + (1/4000)*pos(1)^2 + (1/4000)*pos(2)^2-cos(pos(1)) * cos(sqrt(2)*pos(2)/2);
% end

% Griewank 3d
function fitness = fit_fun(pos, n)
    fitness = 1 + (1/4000)*pos(1)^2 + (1/4000)*pos(2)^2 + (1/4000)*pos(3)^2 - cos(pos(1)) * cos(sqrt(2)*pos(2)/2) * cos(3^(1/3) * pos(3) / 3);
end
