function z = FitnessFunction(x, nRow)
    matrix = reshape(x,[nRow nRow]);
    goldNumber = nRow * ((nRow^2)+1)/2;
    z = 0;
    for i=1:nRow
        z = z + abs( sum(matrix(i,:))-goldNumber );
        z = z + abs( sum(matrix(:,i))-goldNumber );
    end
    z = z + abs( sum(diag(matrix))-goldNumber );
    newMatrix = fliplr(matrix);
    z = z + abs( sum(diag(newMatrix))-goldNumber );
end