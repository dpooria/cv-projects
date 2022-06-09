function parents = ParentSelection(population, numberOfParents)
    individual.Position = [];
    individual.Cost = [];
    parents = repmat(individual,numberOfParents, 1);
    
    for i=1:numberOfParents
        randomIndices = randperm(numel(population),2);
        selectedIndividuals = population(randomIndices);
        [~,bestIndividuals] = sort([selectedIndividuals.Cost]);
        bestParentIndex = bestIndividuals(1);
        parents(i) = selectedIndividuals(bestParentIndex);
    end
end