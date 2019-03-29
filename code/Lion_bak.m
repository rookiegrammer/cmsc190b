classdef Lion < handle
    
    properties
        position=[]
        pbest
        pbestval
        sex='m'
        role=''
    end
    
    methods
        function init(me,position,fit_fun)
            me.position = position;
            me.pbest = position;
            me.pbestval = fit_fun.eval(position);
        end
        
        function offsprings = mate(me,males,muprob,minval,maxval)
            if me.sex ~= 'f'
                offsprings = [];
                return
            end
            
            beta = normrnd(0.5,0.1); % crossover point
            o1 = Lion;
            o2 = Lion;
            
            for i=1:length(me.position) % for all dimension
                mgene = 0;
                for j=1:length(males)
                    mgene = mgene + males(j).position(i);
                end
                mgene = mgene/length(males);
                
                mbeta = 1-beta;
                omutate = 0; % no mutate
                
                if rand(1) <= muprob
                    randvalue = (maxval-minval)*rand()+minval;
                    if rand(1) <= 0.5
                        omutate = 1;
                        o1.position(i,1) = randvalue;
                    else
                        omutate = 2;
                        o2.position(i,1) = randvalue;
                    end
                end
                
                if omutate ~= 1
                    o1.position(i,1) = beta.*me.position(i,1) + mbeta.*mgene;
                end
                
                if omutate ~= 2
                    o2.position(i,1) = mbeta.*me.position(i,1) + beta.*mgene;
                end
            end
            
            % select gender, arrange by sex
            if rand(1) <= 0.5
                o1.sex = 'm';
                o2.sex = 'f';
                offsprings = [o1 o2];
            else
                o2.sex = 'm';
                o1.sex = 'f';
                offsprings = [o2 o1];
            end
            
            
        end
        
        function go_toward(me,other)
            if me.position == other.position
                return % no change
            end
            
            mypos = me.position;
            
            random = rand();
            direction = (other.pbest - mypos); % direction vector times D distance between best and source
            
            mypos = mypos + 2 * random * direction; % push lion to direction, part 1

            normals = null(direction(:).').*norm(direction); % get all vectors normal to direction vector
            deviate = random * pi / 6; % extent to which the point can deviate from the direction
            for z = 1:size(normals, 2) % for all 'dimensions' or direction the point can deviate from, part 2
                mypos = mypos + normals(:,z) .* deviate * (rand() * 2 - 1); % deviate point in direction perpendicular to original direction
                %             ^ this norm vct ^ extnt of deviate ^ -1 to 1 times of deviation
            end
            
            me.position = mypos;
        end
        
        function nmd_mutate(me, lfitval,minval,maxval)
            improve = (me.pbestval - lfitval)/lfitval;
            if rand() < (0.1+min(0.5,improve))
                dlen = length(me.position);
                me.position = minval+(maxval-minval).*rand(dlen, 1);
            end
        end
        
        function did_update = evaluate(me, fit_fun)
            did_update = false;
            if me.position == me.pbest
                return
            end
            
            newbestval = fit_fun.eval(me.position);
            
            if newbestval < me.pbestval
                me.pbest = me.position;
                me.pbestval = newbestval;
                did_update = true;
            end
        end
        
        function print(me, style)
            if me.sex == 'm'
                color = 'b';
            else
                color = 'r';
            end
            dispos(me.pbest, [style color]);
        end
    end
end

