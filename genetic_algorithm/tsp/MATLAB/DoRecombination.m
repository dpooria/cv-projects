function offsprings = DoRecombination(parents, recombinationRate, CostFunction)
    individual.Position = [];
    individual.Cost = [];
    offsprings = repmat(individual,numel(parents),1);
    numberOfGenes = numel(parents(1).Position);

    
    for i=1:2:numel(parents)
       
        randomNumber = randi(100);
        if randomNumber < recombinationRate
              offSpring1 = nan(1,numberOfGenes);
              offSpring2 = nan(1,numberOfGenes);
%             Recombine Parents(i,i+1)
              randomPosition = randi(numberOfGenes);
              offSpring1(1:randomPosition) = parents(i+0).Position(1:randomPosition);
              offSpring2(1:randomPosition) = parents(i+1).Position(1:randomPosition);
              answer1 = ismember(parents(i+1).Position,offSpring1);
              answer2 = ismember(parents(i+0).Position,offSpring2);
              remindedGenes1 = parents(i+0).Position(~answer2);
              remindedGenes2 = parents(i+1).Position(~answer1);
              offSpring1(randomPosition+1:end) = remindedGenes2;
              offSpring2(randomPosition+1:end) = remindedGenes1;
              offsprings(i+0).Position = offSpring1;
              offsprings(i+0).Cost = CostFunction(offSpring1, sqrt(numberOfGenes));
              offsprings(i+1).Position = offSpring2;
              offsprings(i+1).Cost = CostFunction(offSpring2, sqrt(numberOfGenes));
        else
            offsprings(i:i+1) = parents(i:i+1);
        end
        
    end
    
end