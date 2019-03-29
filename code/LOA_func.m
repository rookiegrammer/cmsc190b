function [gbest,gbestval,fitcount] = LOA_func(fit_fun,dimension,population,iterations,space_min,space_max,generated,title,varargin)

if ischar(title)
    cur_date = title;
else
    cur_date = datetime;
    cur_date.Format = 'uuuu-MM-dd-HH-mm-ss';
    cur_date = char(cur_date);
end

% -----------------------
% Main Variables
% -----------------------

print_lions = false;
view_fit_range=[0 6];
graph_function = false;
view_3d_angle = [-45 -45 90];

print_graphics = false;
print_statistics = true;
print_all_fitness = true;
print_end_population = true;

prides_length = 4;

percent_nomad = 0.2;
percent_roam = 0.2;
percent_sex = 0.8;

mating_rate = 0.3;
mutation_prob = 0.2;

immigration_rate = 0.4;

limit = 150000;

nomad_group = Group;
pride_groups = Group.empty(0, prides_length);

% -----------------------

fits = [];
adapt_fun = Fitness_Function;
adapt_fun.function_handle = fit_fun;
adapt_fun.function_argin = varargin;

% -----------------------
% Initialize Lions
% -----------------------

% Generate Initial Population
if isempty(generated)
rand_lio = rand(dimension, population).*(space_max-space_min)+space_min;
else
rand_lio = generated;
end

rand_pop = Lion.empty(0,population);
for i=1:population
    rand_cli = Lion;
    rand_cli.init(rand_lio(:,i), adapt_fun);
    rand_pop(i) = rand_cli;
end

nomd_siz = round(percent_nomad * population);
nomd_ind = randperm(population, nomd_siz);

nomad_group.type = 'n';
nomad_group.maxsize = nomd_siz;

nomd_tmg = Lion.empty(0,nomd_siz);
for i=1:nomd_siz
    nomd_tmg(i) = rand_pop(nomd_ind(i));
end

rand_pop(nomd_ind) = [];

nomd_fms = round((1-percent_sex)*nomd_siz);
nomd_fid = randperm(nomd_siz, nomd_fms);

nomad_group.females = Lion.empty(0, nomd_fms);
for i=1:nomd_fms
    nomd_ths = nomd_tmg(nomd_fid(i));
    nomd_ths.sex = 'f';
    nomad_group.females(i) = nomd_ths;
end

nomd_tmg(nomd_fid) = [];

nomad_group.males = nomd_tmg;
nomad_group.configure();

prid_rem = population-nomd_siz;
prid_szs = ceil(prid_rem/prides_length);
prid_end = prides_length;
for i=1:prid_end
    if i ~= prid_end
        prid_grp = Lion.empty(0,prid_szs);
    else
        prid_grp = Lion.empty(0,length(rand_pop));
    end
    for j=1:prid_szs
        if length(rand_pop) <= 0
            break;
        end
        prid_ths = rand_pop(1);
        prid_grp(j) = prid_ths;
        rand_pop(1) = [];
    end
    prid_gln = length(prid_grp);
    prid_fml = round(prid_gln*percent_sex);
    prid_fin = randperm(prid_gln, prid_fml);

    prid_fgr = Lion.empty(0, prid_fml);
    for j=1:prid_fml
        prid_tfe = prid_grp(prid_fin(j));
        prid_tfe.sex = 'f';
        prid_fgr(j) = prid_tfe;
    end
    prid_grp(prid_fin) = [];
    prid_thg = Group;
    prid_thg.type = 'p';
    prid_thg.females = prid_fgr;
    prid_thg.males = prid_grp;
	prid_thg.maxsize = length(prid_thg.all_lions());
    prid_thg.configure();
    pride_groups(i) = prid_thg;
end

global_best_fitness = nomad_group.lbestval;
global_best = nomad_group.lbest;

for j=1:prides_length
    iter_gpr = pride_groups(j);
    if iter_gpr.lbestval < global_best_fitness
        global_best_fitness = iter_gpr.lbestval;
        global_best = iter_gpr.lbest;
    end
end

if graph_function
    warning('off', 'MATLAB:fplot:NotVectorized')

    if dimension == 1
        base_fig = fplot(@(x) fit_fun(x, 0), [space_min space_max], '-g');
    else
        if dimension == 2
            base_fig = fsurf(@(x,y) fit_fun([x y], 0),[space_min space_max space_min space_max]);
            set(base_fig,'AdaptiveMeshDensity',0,'MeshDensity',60);
        else
            graph_function = false;
        end
    end
    
    if ~print_lions && graph_function
        if print_graphics
            print(['out/loa-graph.png'], '-dpng');
        end
    end
end

if print_lions

    if dimension == 1
        figure(figure)
        axis([space_min space_max view_fit_range]);
        xlabel('x_1')
        ylabel('f(X)')
    else
        if dimension == 2
            figure(figure)
            axis([space_min space_max space_min space_max view_fit_range]);
            view(view_3d_angle); % 3D angle
            xlabel('x_1')
            ylabel('x_2')
            zlabel('f(X)')
        end
    end
end

if print_lions

    hold on;

    % PRINT ALL
    for j=1:prides_length
        iter_gpr = pride_groups(j);
        iter_gpr.print();
    end
    nomad_group.print();
    fprintf('- Best Fitness: %g\n', global_best_fitness);

    if graph_function
        copyobj(base_fig,gca);
    end

    if print_graphics
        print(['out/loa-iter-0.png'], '-dpng');
    end

    hold off;
end

if print_statistics
    file_id = fopen(['out/loa-iter-stats-' cur_date '.txt'], 'wt');
    fprintf(file_id, 'Iterations: %g\n', iterations);
    fprintf(file_id, 'Function Evaluations Limit: %g\n', limit);
    fprintf(file_id, 'Population: %g\n\n', population);
    fprintf(file_id, 'Prides Length: %g\n\n', prides_length);
    fprintf(file_id, 'Percent Nomads: %g\n', percent_nomad);
    fprintf(file_id, 'Percent Roaming: %g\n', percent_roam);
    fprintf(file_id, 'Percent Sex: %g\n\n', percent_sex);
    fprintf(file_id, 'Mating Rate: %g\n', mating_rate);
    fprintf(file_id, 'Mutation Probability: %g\n', mutation_prob);
    fprintf(file_id, 'Immigration Rate: %g\n\n', immigration_rate);
    fprintf(file_id, 'Dimensions: %g\n', dimension);
    fprintf(file_id, 'N-Dimension Space: %g:%g\n', space_min, space_max);
    fprintf(file_id, '\nFitness Iterations:\n%g', global_best_fitness);
end

% -----------------------
% Start Generations
% -----------------------

if print_statistics
        if print_all_fitness
          fprintf(file_id, '\n');
          lfits = [];
          for j=1:prides_length
            pfits = pride_groups(j).fitnesses();
            fprintf(file_id, '%g, ', pfits );
            lfits = [lfits pfits];
          end
          nfits = nomad_group.fitnesses();
          fprintf(file_id, '%g, ', nfits );
          lfits = [lfits nfits];
          
          fits = [fits; lfits NaN(1,size(fits,2)-length(lfits))];
          fprintf(file_id, '\n' );
        end
end

for i=1:iterations

    if print_lions
        fprintf('Iteration %d\n', i);
    end

    for j=1:prides_length
        iter_gpr = pride_groups(j);
        iter_gpr.recount();
        iter_gpr.do_pride_fem(percent_roam,space_min,space_max, adapt_fun);
        iter_gpr.do_pride_mal(percent_roam, adapt_fun,space_min,space_max);
        iter_gpr.mate(mating_rate,mutation_prob,space_min,space_max,adapt_fun);
        iter_gpr.equilibriate(nomad_group,percent_sex);
        pride_groups(j) = iter_gpr;
    end

    nomad_group.recount();
    nomad_group.do_nomad_all(space_min, space_max, adapt_fun);
	nomad_group.mate(mating_rate,mutation_prob,space_min,space_max,adapt_fun);
	nomad_group.invade(pride_groups);

    for j=1:prides_length
		iter_gpr = pride_groups(j);
		iter_gpr.emigrate(percent_sex,immigration_rate,nomad_group);
    end

 	nomad_group.immigrate(pride_groups,percent_sex);
	nomad_group.equilibriate(nomad_group,percent_sex);

    % CHECK FITNESS
    if nomad_group.lbestval < global_best_fitness
        global_best_fitness = nomad_group.lbestval;
        global_best = nomad_group.lbest;
    end
    for j=1:prides_length
        iter_gpr = pride_groups(j);
        if iter_gpr.lbestval < global_best_fitness
            global_best_fitness = iter_gpr.lbestval;
            global_best = iter_gpr.lbest;
        end
    end

    if print_statistics
        if print_all_fitness
          fprintf(file_id, '\n');
          lfits = [];
          for j=1:prides_length
            pfits = pride_groups(j).fitnesses();
            fprintf(file_id, '%g, ', pfits );
            lfits = [lfits pfits];
          end
          nfits = nomad_group.fitnesses();
          fprintf(file_id, '%g, ', nfits );
          lfits = [lfits nfits];
          
          fits = [fits; lfits NaN(1,size(fits,2)-length(lfits))];
          fprintf(file_id, '\n' );
        else
          fprintf(file_id, '\n%g', global_best_fitness);
        end
    end

    if print_lions

        cla
        hold on;

        % PRINT ALL
        for j=1:prides_length
            iter_gpr = pride_groups(j);
            iter_gpr.print();
        end
        nomad_group.print();
        fprintf('- Best Fitness: %g\n', global_best_fitness);

        if graph_function
            copyobj(base_fig,gca);
        end

        if print_graphics
            print(['out/loa-iter-' num2str(i) '.png'], '-dpng');
        end

        hold off;
        pause(0.0001)
    end

    if adapt_fun.fitness_count >= limit
        break;
    end

end

if print_statistics
    fprintf(file_id, '\nBest:\n');
    fprintf(file_id, '%g\t', global_best);
    fprintf(file_id, '\n');
    fprintf(file_id, '\nFitness Evaluations:\n%g\n', adapt_fun.fitness_count);
    
    if print_end_population
      fprintf(file_id, '\nEnd Population\n');
      for i=1:prides_length
        fprintf(file_id, '(%s), ', [pride_groups(i).all_lions().pbest] );
      end
      fprintf(file_id, '(%s), ', [nomad_group.all_lions().pbest] );
    end
    fclose(file_id);
end

gbest = global_best';
gbestval = global_best_fitness;
fitcount = adapt_fun.fitness_count;

csvwrite(['out/loa-fits-' cur_date '.csv'], fits);

end
