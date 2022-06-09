function offsprings = DoMutation(offsprings, mutationRate, CostFunction)
    numberOfGenes = numel(offsprings(1).Position);

    for i=1:numel(offsprings)
        randomNumber = randi(100,1,1);
        if randomNumber < mutationRate
            % Do Mutation
            randomPositions = zeros(2,1);
            while randomPositions(1) == randomPositions(2)
                randomPositions = randi(numberOfGenes,2,1);
            end
            offsprings(i).Position([randomPositions(1) randomPositions(2)]) = offsprings(i).Position([randomPositions(2) randomPositions(1)]);
            offsprings(i).Cost = CostFunction(offsprings(i).Position, sqrt(numberOfGenes));
        end
    end
end