clc;
clear;
close all;

%% Define Problem
numberofCities = 25;

%% Parameters
numberOfPopulation = 300;
numberOfGenerations = 500;
numberOfParents = floor(numberOfPopulation / 2);
recombinationRate = 100;
mutationRate = 30;

%% Initialization
City.label = [];
City.x = [];
City.y = [];
distanceOfCities = @(c1, c2) (sqrt((c1.x - c2.x).^2 + (c1.y - c2.y).^2));

Cities = repmat(City, numberOfCities, 1);
individual.Position = nan(1,numberOfGenes);
individual.Cost = nan(1,1);
population = repmat(individual,numberOfPopulation, 1);
parents = repmat(individual,numberOfParents, 1);
offsprings = repmat(individual,numberOfParents, 1);
bestIndividuals = nan(numberOfGenerations + 1, 1);
for i=1:numberOfPopulation
    population(i).Position = randperm(numberOfGenes);
    population(i).Cost = CostFunction(population(i).Position,numberOfRows);
end

%% Main Loop
PopulationBeforeEvolution = population;
[~, sortOrder] = sort([PopulationBeforeEvolution.Cost]);
PopulationBeforeEvolution = PopulationBeforeEvolution(sortOrder);
bestIndividuals(1) = PopulationBeforeEvolution(1).Cost;
tic;
for i=1:numberOfGenerations
    parents = ParentSelection(population, numberOfParents);
    offsprings = DoRecombination(parents, recombinationRate, CostFunction);
    offsprings = DoMutation(offsprings, mutationRate, CostFunction);
    [population, bestIndividuals(i+1)] = DoSurvivorSelection(population, offsprings);
    disp(['Generatios #' num2str(i), '  Best Individual Fitness: ' num2str(bestIndividuals(i+1))]);
end
elapsedTime = toc;
disp('**********************************************');
disp(['Elapsed Time = ' num2str(elapsedTime)]);
disp(['Best Individual = ' mat2str(reshape(population(1).Position,[numberOfRows numberOfRows])) ' ,Fitness = ' num2str(bestIndividuals(end))]);

%% Results
figure;
subplot(1,2,1);
plot([PopulationBeforeEvolution.Cost]);
title('Populations Cost Before Evolution');
subplot(1,2,2);
plot([population.Cost]);
title('Populations Cost After Evolution');

figure;
plot(bestIndividuals);
title('Fitness-Generation');
xlabel('Generations');
ylabel('Fitness');
