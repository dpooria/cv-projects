function [population, bestIndividual] = DoSurvivorSelection(population, offsprings)
    numberOFPopulation = numel(population);
    newPopulation = [population;offsprings];
    [~, sortOrder] = sort([newPopulation.Cost]);
    newPopulation = newPopulation(sortOrder);
    population = newPopulation(1:numberOFPopulation);
    bestIndividual = population(1).Cost;
end