classdef iLOA_setting < handle
    properties
        numberOfPrides
        percentNomad
        percentRoam
        percentSex
        rateMating
        probabilityMutation
        rateImmigration
        percentGroupInfluence
        annealing
        pressureRankedSelect
        pressureNearBest
    end
    methods
        function params = get_all(me)
            params = [me.numberOfPrides me.percentNomad me.percentRoam me.percentSex me.rateMating me.probabilityMutation me.rateImmigration me.percentGroupInfluence me.annealing me.pressureRankedSelect me.pressureNearBest];
        end
    end
end

