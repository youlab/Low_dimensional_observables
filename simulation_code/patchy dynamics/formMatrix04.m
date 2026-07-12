% Author: Emrah Simsek

function [F, G] = formMatrix04(n1, n2, interactionVectors)

    PosIntVectors = interactionVectors;
    
    NegIntVectors = interactionVectors;
    
    for i = 1:n1
    
        for j = 1:2*n2
    
              M = PosIntVectors{i, j}; N = NegIntVectors{i, j};
              M(M < 0) = 0;  N(N > 0) = 0;  
              PosIntVectors{i, j} = M; NegIntVectors{i, j} = -N;
    
        end
    
    end

    % ix = (col - 1) * n + row; %linear index

    % define cooperation (positive interaction) input matrix
    F = zeros(2*n1*n2, 2*n1*n2); 

    % define competition (negative interaction) input matrix
    G = zeros(2*n1*n2, 2*n1*n2); 

    k = 1;

    % assign the elements of these matrices
    for j = 1:2*n2

        for i = 1:n1

      % kth row of F defines the positive input to the native population at the {i, j} site or focal populations sitting at the (k - n1*n2) site from
        % others
            F(k, :) = PosIntVectors{i, j};

      % kth row of G defines the positive input to the native population at the {i, j} site or focal populations sitting at the (k - n1*n2) site from
        % others
            G(k, :) = NegIntVectors{i, j};

            k = k + 1;

        end

    end
    
end