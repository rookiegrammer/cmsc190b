function [gbest,gbestval,fitcount,iter] = iLOA_func(fit_fun,dimension,population,iterations,space_min,space_max,generated,title,setting,varargin)

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

% ________OPTIONS________
convergence_mode = true; %#ok<*NASGU>
convergence_tolerance = 0.001;

view_fit_range   = [0 5];
view_3d_angle    = [-45 -45 90];

graph_function = false;
print_lions = false;
print_graphics = false; % Save Iterations Images to File
print_statistics = false; % Save Statistics per Iteration to file
print_all_fitness = false;
print_end_population = false;

limit = 150000;

% ______PARAMETERS_______
prides_length = 4;

percent_nomad = 0.2;
percent_roam = 0.8;
percent_sex = 0.2;

mating_rate = 0.6;
mutation_prob = 0.8;

immigration_rate = 0.2;

% ____NEW_PARAMETERS______
percent_influence = 0.8;
do_annealing = true;
selection_pressure = 2;
nearness_pressure = 2;

% ______USERS_END_________
% ________________________

% Load object based passthrough settings
if isobject(setting)
    prides_length = setting.numberOfPrides;
    percent_nomad = setting.percentNomad;
    percent_roam = setting.percentRoam;
    percent_sex = setting.percentSex;
    mating_rate = setting.rateMating;
    mutation_prob = setting.probabilityMutation;
    immigration_rate = setting.rateImmigration;
    percent_influence = setting.percentGroupInfluence;
    do_annealing = setting.annealing;
    selection_pressure = setting.pressureRankedSelect;
    nearness_pressure = setting.pressureNearBest;
end

% -----------------------
% Initialize necessary variables
fits = [];
nomad_group = iGroup;
nomad_group.anneal = do_annealing;
pride_groups = iGroup.empty(0, prides_length);
file_id = '';

% -----------------------

% Create Fitness Function Object
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

% Initialize First Lion Objects
rand_pop = iLion.empty(0,population);
for i=1:population
    rand_cli = iLion;
    rand_cli.init(rand_lio(:,i), adapt_fun);
    rand_pop(i) = rand_cli;
end

% Calculate How Many Nomads and Create Random Indices
nomd_siz = round(percent_nomad * population);
nomd_ind = randperm(population, nomd_siz);

% Set Nomad group type and size limit
nomad_group.type = 'n';
nomad_group.maxsize = nomd_siz;

% Create Lion Array for Nomads
nomd_tmg = iLion.empty(0,nomd_siz);
for i=1:nomd_siz
    nomd_tmg(i) = rand_pop(nomd_ind(i));
end

% Remove from Original Array
rand_pop(nomd_ind) = [];

% Calculate Random Female Nomad Indices
nomd_fms = round((1-percent_sex)*nomd_siz);
nomd_fid = randperm(nomd_siz, nomd_fms);

% Create Array for Female Nomads with initialization
nomad_group.females = iLion.empty(0, nomd_fms);
for i=1:nomd_fms
    nomd_ths = nomd_tmg(nomd_fid(i));
    nomd_ths.sex = 'f';
    nomad_group.females(i) = nomd_ths;
end

% Remove from original Nomad Array
nomd_tmg(nomd_fid) = [];

% All Remaining Nomad Lions are Male
nomad_group.males = nomd_tmg;
nomad_group.configure();

% Get Partition Size for Prides
prid_rem = population-nomd_siz;
prid_szs = ceil(prid_rem/prides_length);
prid_end = prides_length;

% Create Respective Prides
for i=1:prid_end
    % Truncate Remaining to a pride
    if i ~= prid_end
        prid_grp = iLion.empty(0,prid_szs);
    else
        prid_grp = iLion.empty(0,length(rand_pop));
    end

    % Put to pride
    for j=1:prid_szs
        if length(rand_pop) <= 0
            break;
        end
        prid_ths = rand_pop(1);
        prid_grp(j) = prid_ths;
        rand_pop(1) = [];
    end

    % Determine Pride Gender Sizes
    prid_gln = length(prid_grp);
    prid_fml = round(prid_gln*percent_sex);
    prid_fin = randperm(prid_gln, prid_fml);

    % Divide Pride Population to Male and Female to their placeholder
    prid_fgr = iLion.empty(0, prid_fml);
    for j=1:prid_fml
        prid_tfe = prid_grp(prid_fin(j));
        prid_tfe.sex = 'f';
        prid_fgr(j) = prid_tfe;
    end
    prid_grp(prid_fin) = [];
    prid_thg = iGroup;
    prid_thg.type = 'p';
    prid_thg.females = prid_fgr;
    prid_thg.males = prid_grp;

    % Configure Pride
	prid_thg.maxsize = length(prid_thg.all_lions());
    prid_thg.configure();
    pride_groups(i) = prid_thg;
end

% Find best fitness in nomad group
global_best_fitness = nomad_group.lbestval;
global_best = nomad_group.lbest;

% Find best fitness in prides
for j=1:prides_length
    iter_gpr = pride_groups(j);
    if iter_gpr.lbestval < global_best_fitness
        global_best_fitness = iter_gpr.lbestval;
        global_best = iter_gpr.lbest;
    end
end

% Graph function for the first time
if graph_function
    warning('off', 'MATLAB:fplot:NotVectorized')
    
    if dimension > 2
        graph_function = false;
    end
    
    do_graph_function();

    if ~print_lions && graph_function
        if print_graphics
            print(['out/iloa-graph-' title '.png'], '-dpng'); %#ok<*UNRCH>
        end
    end
end

% Make axis for the lions
if print_lions
    do_make_axis();
end

% Do a first time plot or surface for the function and lions in canvas
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
        do_graph_function();
    end

    if print_graphics
        print(['out/iloa-iter-0-' cur_date '.png'], '-dpng');
    end

    hold off;
end

% Start creating info text file
if print_statistics
    LogicalStr = {'false', 'true'};
    file_id = fopen(['out/iloa-iter-stats-' cur_date '.txt'], 'wt');
    fprintf(file_id, 'Iterations: %g\n', iterations);
    fprintf(file_id, 'Function Evaluations Limit: %g\n', limit);
    fprintf(file_id, 'Population: %g\n\n', population);
    fprintf(file_id, 'Prides Length: %g\n\n', prides_length);
    fprintf(file_id, 'Percent Nomads: %g\n', percent_nomad);
    fprintf(file_id, 'Percent Roaming: %g\n', percent_roam);
    fprintf(file_id, 'Percent Sex: %g\n\n', percent_sex);
    fprintf(file_id, 'Mating Rate: %g\n', mating_rate);
    fprintf(file_id, 'Mutation Probability: %g\n', mutation_prob);
    fprintf(file_id, 'Immigration Rate: %g\n', immigration_rate);
    fprintf(file_id, '\nPercent Group Influence: %g\n', percent_influence);
    fprintf(file_id, 'Annealing On: %s\n', LogicalStr{do_annealing+1} );
    fprintf(file_id, 'Ranked Selection Pressure: %g\n', selection_pressure);
    fprintf(file_id, 'Near To Best Random Pressure: %g\n\n', nearness_pressure);
    fprintf(file_id, 'Dimensions: %g\n', dimension);
    fprintf(file_id, 'N-Dimension Space: %g:%g\n', space_min, space_max);
    fprintf(file_id, '\nFitness Iterations:\n%g (Initial Best)\n', global_best_fitness);
end

% -----------------------
% Start Generations
% -----------------------

% Print lions in info text
if print_statistics && print_all_fitness
    do_print_fitness();
end

% Start the loop until max iteration
for i=1:iterations

    if print_lions
        fprintf('Iteration %d\n', i);
    end
    
    % Process all pride lions
    for j=1:prides_length
        iter_gpr = pride_groups(j);
        iter_gpr.recount();
        iter_gpr.do_pride_fem(percent_roam,space_min,space_max, adapt_fun, percent_influence,selection_pressure);
        iter_gpr.do_pride_mal(percent_roam, adapt_fun,space_min,space_max, percent_influence,selection_pressure);
    end
    
    for j=1:prides_length
        iter_gpr = pride_groups(j);
        iter_gpr.mate(mating_rate,mutation_prob,space_min,space_max,adapt_fun,selection_pressure);
        iter_gpr.equilibriate(nomad_group,percent_sex);
    end
    
    % Process nomad group next
    nomad_group.recount();
    nomad_group.do_nomad_all(space_min, space_max, adapt_fun, global_best, nearness_pressure);
  	nomad_group.mate(mating_rate,mutation_prob,space_min,space_max,adapt_fun,selection_pressure);
  	nomad_group.invade(pride_groups);
    
    % Continue with prides emigration
    for j=1:prides_length
		iter_gpr = pride_groups(j);
		iter_gpr.emigrate(percent_sex,immigration_rate,nomad_group);
    end
    
    % Finish with nomad_group 
 	nomad_group.immigrate(pride_groups,percent_sex);
	nomad_group.equilibriate(nomad_group,percent_sex);

    % CHECK FITNESS Get minimum
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

    % Print to info text
    if print_statistics && print_all_fitness
        do_print_fitness();
    end
    
    % Create preview graph
    if print_lions

        cla
        hold on;

        % PRINT ALL
        for j=1:prides_length
            iter_gpr = pride_groups(j);
            iter_gpr.print();
            fprintf('-- Pride %d: %d + %d\n', j, length(iter_gpr.males), length(iter_gpr.females));
        end
        nomad_group.print();
        fprintf('-- Nomads : %d + %d\n', length(nomad_group.males), length(nomad_group.females));

        fprintf('- Best Fitness: %g\n', global_best_fitness);

        if graph_function
            do_graph_function();
        end

        if print_graphics
            print(['out/iloa-iter-' num2str(i) '-' cur_date '.png'], '-dpng');
        end

        hold off;
        pause(0.0001)
    end
    
    if check_converge()
        break;
    end
    
    % Check limit evaluation count
    if adapt_fun.fitness_count >= limit
        break;
    end
    
    

end

% Append to text with summary
if print_statistics
    fprintf(file_id, '\nBest:\n(');
    fprintf(file_id, '%g, ', global_best);
    fprintf(file_id, ')\n');
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

% Generate output
gbest = global_best';
gbestval = global_best_fitness;
fitcount = adapt_fun.fitness_count;
iter = i;

if ~isempty(fits)
    csvwrite(['out/iloa-fits-' cur_date '.csv'], fits);
end


    function converged = check_converge()
        pbestvals = [];
        for z=1:prides_length
            pbestvals = [pbestvals [pride_groups(z).all_lions().pbestval]];
        end
        vals = abs(max(pbestvals)-min(pbestvals));
        converged = vals <= convergence_tolerance;
        
    end

    function do_graph_function()
        if dimension == 1
            fplot(@(x) fit_fun(x, 0), [space_min space_max], '-g');
        else
            if dimension == 2
                base_fig = fsurf(@(x,y) fit_fun([x y], 0),[space_min space_max space_min space_max]);
                set(base_fig,'AdaptiveMeshDensity',0,'MeshDensity',60);
            end
        end
    end

    function do_make_axis()
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

    function do_print_fitness()
        fprintf(file_id, '\n');
        lfits = [];
        for z=1:prides_length
            pfits = pride_groups(z).fitnesses();
            fprintf(file_id, '%g, ', pfits );
            lfits = [lfits pfits]; %#ok<*AGROW>
        end
        nfits = nomad_group.fitnesses();
        fprintf(file_id, '%g, ', nfits );
        lfits = [lfits nfits];

        fits = [fits; lfits NaN(1,size(fits,2)-length(lfits))];
        fprintf(file_id, '\n' );
    end

end


