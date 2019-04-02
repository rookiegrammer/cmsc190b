classdef iGroup < handle
    properties
        type

        lbest % local best
        lbestval % local best value

        lworstval
        lcurworstval = -Inf;

        anneal = false;

        males = [];
        females = [];

        improved = 0;

		maxsize = 0;
    end
    methods
        function lions = all_lions(me)
            lions = [me.males me.females];
        end

        function recount(me)
            me.improved = 0;
            if me.anneal
                me.lworstval = me.lcurworstval;
                me.lcurworstval = -Inf;
            end
        end

        function configure(me)
            lions = me.all_lions();
            me.lbestval = Inf;
            for i=1:length(lions)
                lion_ith = lions(i);
                if lion_ith.pbestval < me.lbestval
                    me.lbest = lion_ith.pbest;
                    me.lbestval = lion_ith.pbestval;
                end
                if me.anneal && lion_ith.pbestval > me.lcurworstval
                    me.lcurworstval = lion_ith.pbestval;
                end
            end
            if me.anneal
                me.lworstval = me.lcurworstval;
            end
        end

        function force_get_best(me)
            lions = me.all_lions();
            pbestval = me.lbestval;
            for i=1:length(lions)
                if lions(i).pbestval < pbestval
                    pbest = lions(i).pbest;
                    pbestval = lions(i).pbestval;
                end
            end
            me.lbest = pbest;
            me.lbestval = pbestval;
        end

        function do_pride_fem(me, roam_per, min_val, max_val, fit_fun, grp_inf, sel_press)
                fem_len = length(me.females);

                if fem_len == 0
                    return
                end

                roam_len = floor(roam_per * fem_len);
                roam_ind = randperm(fem_len, roam_len);

                hunt_dim = length(me.females(1).position);
                hunt_ers = me.females;
                hunt_ers( roam_ind ) = [];

                hunt_len = length(hunt_ers);
                hunt_cln = ceil(hunt_len / 3);

                [~, indices]=sort([hunt_ers.pbestval]);
                for i=1:hunt_cln
                    hunt_ers(indices(i)).role = 'c';
                end

                hunt_pry = min_val + rand(hunt_dim,1)*(max_val-min_val);

                for i=1:hunt_len
                    hunt_cli = hunt_ers(i);
                    hunt_rnd = rand();
                    hunt_lst = hunt_cli.pbestval;
                    hunt_lps = hunt_cli.position;

                    if hunt_ers(i).role == 'c'
                        hunt_cli.position = hunt_lps + hunt_rnd .* (hunt_pry - hunt_lps);
                    else
                        hunt_cli.position = hunt_pry + hunt_rnd .* (hunt_pry - hunt_lps);
                    end

                    hunt_cli.evaluate(fit_fun, min_val, max_val);

                    if hunt_cli.pbestval < hunt_lst
                        hunt_pim = (hunt_lst - hunt_cli.pbestval) / hunt_lst;
                        hunt_prn = rand();
                        hunt_pry = hunt_pry + hunt_prn .* hunt_pim .* (hunt_pry - hunt_cli.position);
                        me.improved = me.improved + 1;

                        me.update_best(hunt_cli);
                    end
                end

                roam_all = me.all_lions();
                [~, roam_srt] = sort([roam_all.pbestval]);
                roam_all = roam_all(roam_srt);

                roam_tsiz = max(2, ceil(me.improved/2));

                for i=1:length(roam_ind)
                    roam_cli = me.females(roam_ind(i));
                    roam_tin = randk(length(roam_all), roam_tsiz, sel_press);
                    roam_tim = false;
                    for j=1:roam_tsiz
                        roam_oth = roam_all(roam_tin(j));

                        roam_clf = roam_cli.pbestval;

                        roam_cli.go_toward(roam_oth, me.lbest, grp_inf);
                        roam_cli.evaluate(fit_fun, min_val, max_val);

                        if roam_cli.pbestval < roam_clf
                            roam_tim = true;
                            me.update_best(roam_cli);
                        end
                    end
                    roam_cli.position = roam_cli.pbest;
                    if roam_tim
                        me.improved = me.improved + 1;
                    end
                end
        end

        function do_pride_mal(me, roam_per, fit_fun, min_val, max_val, grp_inf, sel_press)
            male_len = length(me.males);

            roam_all = me.all_lions();
            [~, roam_srt] = sort([roam_all.pbestval]);
            roam_all = roam_all(roam_srt);

            roam_len = ceil(roam_per * male_len);

            if roam_len > 0
                for i=1:male_len
                    roam_cli = me.males(i);
                    roam_tin = randk(length(roam_all), roam_len, sel_press);
                    roam_tim = false;
                    for j=1:roam_len
                        roam_oth = roam_all(roam_tin(j));

                        roam_clf = roam_cli.pbestval;

                        roam_cli.go_toward(roam_oth, me.lbest, grp_inf);
                        roam_cli.evaluate(fit_fun, min_val, max_val);

                        if roam_cli.pbestval < roam_clf
                            roam_tim = true;
                            me.update_best(roam_cli);
                        end
                    end
                    roam_cli.position = roam_cli.pbest;
                    if roam_tim
                        me.improved = me.improved + 1;
                    end
                end
            end
        end

        function do_nomad_all(me, min_val, max_val, fit_fun, g_best, near_press)
            nomd_all = me.all_lions();
            nomd_len = length(nomd_all);
            nomd_bft = me.lbestval;
            for i=1:nomd_len
                nomd_cli = nomd_all(i);
                nomd_lft = nomd_cli.pbestval;
                nomd_imp = (nomd_lft - nomd_bft) / nomd_bft;
                nomd_prb = 0.1 + min(0.5, nomd_imp);

                if rand() <= nomd_prb
                    nomd_cli.random_value(g_best, min_val, max_val, near_press);
                    nomd_cli.evaluate(fit_fun, min_val, max_val);
                    if nomd_cli.pbestval < nomd_lft
                        me.improved = me.improved + 1;
                        me.update_best(nomd_cli);
                    end
                end
            end
        end

		function invade(me,pride_grps)
			if me.type == 'p'
				return;
			end
			%nomad defense
			for i=1:length(me.males) % for each nomads
				for j=1:length(pride_grps)
					if rand(1)<=0.5 % 50% probability that nomad will attack pride
						nm = me.males(i); % nomad male
						ipride = pride_grps(j);% invaded pride
                        success = false;
						for k=1:length(ipride.males)
							rm = ipride.males(k); % resident male
							if rm.pbestval > nm.pbestval % true if nomad male is stronger than resident
								% switch places of resident and nomad
								ipride.males(k) = nm;
                                ipride.update_best(nm);

								me.males(i) = rm;
                                me.update_best(rm);

                                success = true;
								break;
							end
                        end
                        if success
                            break; % the nomad has invaded one pride now we escape
                        end
					end
				end
			end
        end

        % PRIDE -> NOMAD
		function emigrate(me,sex_rate,im_rate,nomad_grp)
			if me.type == 'n'
				return;
            end
            
            femlen = length(me.females);
            
            if femlen == 0
                return
            end

			maxfem = sex_rate*me.maxsize; %maximum females in this pride
			imfem = fix(im_rate*maxfem); % immigrating females

			mifem = fix(femlen-maxfem + imfem);% number of migrating females (surplus + migrating)
            
            if mifem == femlen
                mifem = mifem - 1;% stop algorithm from disposing all females
            end
            
			imifem = randperm(length(me.females),mifem); %indices of migrating females

            nomad_grp.add_new(me.females(imifem)); % get migrating...

            me.females(imifem)=[]; % remove them here
		end

		function mate(me,mating_rate,mutation_prob,space_min,space_max,fit_fun, sel_press)
            tgrp_mln = length(me.males);
            tgrp_fln = length(me.females);
            
            if tgrp_mln == 0 || tgrp_fln == 0
                return
            end

            [~, ind] = sort([me.females.pbestval]);
            me.females = me.females(ind);

            [~, ind] = sort([me.males.pbestval]);
            me.males = me.males(ind);

			fheat = fix(mating_rate * tgrp_fln);% number of females in heat
			ifheat = randk(tgrp_fln, fheat, sel_press); % indices of females in heat
            
            tgrp_typ = me.type;
            tgrp_mht = randi(tgrp_mln); % index/number of male/s in heat

            chld_m = iLion.empty(0, fheat); % prevent offsprings get r4ped
            chld_f = iLion.empty(0, fheat);

            for i=1:fheat
                tgrp_fem = me.females(ifheat(i));
                if tgrp_typ == 'p'
                    imheat = randk(tgrp_mln,tgrp_mht, sel_press); % indices of males in heat
                    mheat = me.males(imheat);
                    offsprings = tgrp_fem.mate(mheat, mutation_prob, space_min, space_max); % mate
                else
                    tgrp_mht = randk(tgrp_mln, 1, sel_press);
                    offsprings = tgrp_fem.mate(me.males(tgrp_mht), mutation_prob, space_min, space_max);
                end

                %set offspring fitness and compare
				offsprings(1).init(offsprings(1).position,fit_fun);
                offsprings(2).init(offsprings(2).position,fit_fun);

                % index 1 is standard for the male
				chld_m(i) = offsprings(1);
				chld_f(i) = offsprings(2);
            end

            me.add_new(chld_m);
            me.add_new(chld_f);
        end

        function equilibriate(me,nomad_group,sex_rate)
            if me.type == 'p'
                %internal defense
				todrout = (length(me.males))-ceil((1-sex_rate)*me.maxsize);%number of males to drive out (internal defense)

				%sort males from weakest to strongest
				[~, ind] = sort([me.males.pbestval],'descend');
                out_indices = ind(1:todrout); % get first n indices of weakest
                nomad_group.add_new(me.males(out_indices)); % add to nomad males with indices
                me.males(out_indices) = []; % remove them here
            else
                %kill the weak
				mtokill = floor(length(me.males)-me.maxsize*sex_rate); % number of males to kill
				ftokill = floor(length(me.females)-me.maxsize*(1-sex_rate));% number of females to kill

				%sort nomad females from weakest to strongest
				[~, ind] = sort([me.females.pbestval],'descend');
                weak_indices = ind(1:ftokill);
				me.females(weak_indices) = [];

				%sort nomad males from weakest to strongest
				[~, ind] = sort([me.males.pbestval],'descend');
                weak_indices = ind(1:mtokill);
				me.males(weak_indices) = [];
            end
        end

        % NOMAD -> EMPTY Prides
		function immigrate(me,pride_grps,sex_rate)
			if me.type == 'p'
				return;
            end

            %sort females from strongest to get the strong ones first
            [~, ind] = sort([me.females.pbestval]);
            me.females = me.females(ind);

            prid_len = length(pride_grps);

            need_fem = zeros(1, prid_len);

            % how much do each pride need?
            for i=1:prid_len
				pgrp = pride_grps(i);
                need_fem(i) = max(0, fix(sex_rate * pgrp.maxsize + 0.0001) - length(pgrp.females)); % do you need some?
            end

            need_len = sum(need_fem);
            if need_len > length(me.females)
                return; % Don't take anymore if need is greater
            end

            ifem = randperm(need_len); % shuffle indices 1 to needed
            j = 0;

			% give each pride the strong ones
			for i=1:prid_len
                pgrp = pride_grps(i);

                many = need_fem(i); % females this pride needs
                k = j + many;
                first_ind = ifem(j+1:k); % get next n indices

                pgrp.add_new(me.females(first_ind)); % give it to the pride

                j = k;
            end
            me.females(ifem)=[]; % clear those at my indices
        end

        function update_best(me, lion)
            if lion.pbestval < me.lbestval
                me.lbestval = lion.pbestval;
                me.lbest = lion.pbest;
            else
                if me.anneal && lion.pbestval > me.lcurworstval
                    me.lcurworstval = lion.pbestval;
                end
            end
        end

        function add_new(me, lions)
            for i=1:length(lions)
                lion = lions(i);
                add = true;
                if me.anneal && lion.pbestval > me.lworstval
                    add = false;
                end
                if add
                    if lion.sex == 'm'
                        me.males = [me.males lion];
                    else
                        me.females = [me.females lion];
                    end
                    me.update_best(lion);
                end
            end
        end

        function array = fitnesses(me)
          lions = me.all_lions();
          array = [lions.pbestval];
        end

        function print(me)
%             return; %disable print

            if me.type == 'p'
                style = '*';
            else
                style = '+';
            end

            len_m = length(me.males);
            len_f = length(me.females);

            for i=1:len_m
                me.males(i).print(style);
            end
            for i=1:len_f
                me.females(i).print(style);
            end
        end
    end
end
