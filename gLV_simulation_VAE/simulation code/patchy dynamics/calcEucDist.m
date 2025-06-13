function gridDistances = calcEucDist(n1, n2, dmin)

    % Create a matrix to store pairwise Euclidean distances
    gridDistances = zeros(n1*n2, n1*n2);
    
    % Compute pairwise Euclidean distances
    for i = 1:n1*n2
        for j = 1:n1*n2
            [row_i, col_i] = ind2sub([n1, n2], i);
            [row_j, col_j] = ind2sub([n1, n2], j);
            
            distance = norm([row_i, col_i] - [row_j, col_j]);
            gridDistances(i, j) = dmin * distance;
        end
    end

end