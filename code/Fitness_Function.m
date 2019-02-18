classdef Fitness_Function < handle
    
    properties
        function_argin
        function_handle
        fitness_count = 0;
    end
    
    methods
        function output = eval(me, vector)
            me.fitness_count = me.fitness_count + 1;
            output = feval(me.function_handle, vector, me.function_argin{:});
        end
    end
end

