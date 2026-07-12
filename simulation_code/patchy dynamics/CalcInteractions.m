function [interactionVectors, Nat_growthRateMatrix, Nf_growthRateMatrix, gridDistances] = CalcInteractions(n1, n2, gridDistances, Nat_speciesAssignments, Nf_speciesAssignments, interactionMatrix, muVec)

    % Initialize interaction vectors
    interactionVectors = cell(n1, 2 * n2);
    
    Nat_growthRateMatrix =zeros(n1, n2);
    Nf_growthRateMatrix = zeros(n1, n2);

    % Compute interaction vectors and growth rate matrix for each grid site
    for row = 1:n1
    
        for col = 1:n2
    
            % Scale interactions based on distances
            distance_indices = sub2ind([n1, n2], row, col);
            
            gridDistances(gridDistances == 0) = 1; % maximize the local populations' self-interactions
            
            Nat_species = Nat_speciesAssignments(row, col);
            Nf_species = Nf_speciesAssignments(row, col);
    
            % Get interaction strengths from the interaction matrix for the selected species
            NatNat_interactions = interactionMatrix(Nat_species, Nat_speciesAssignments);
            NfNat_interactions = interactionMatrix(Nat_species, Nf_speciesAssignments);
            NatNf_interactions = interactionMatrix(Nf_species, Nat_speciesAssignments);
            NfNf_interactions = interactionMatrix(Nf_species, Nf_speciesAssignments);
    
            scaled_NatNat_Interactions = NatNat_interactions ./ gridDistances(distance_indices, :);
            scaled_NfNat_Interactions = NfNat_interactions ./ gridDistances(distance_indices, :);
            scaled_NatNf_Interactions = NatNf_interactions ./ gridDistances(distance_indices, :);
            scaled_NfNf_Interactions = NfNf_interactions ./ gridDistances(distance_indices, :);
    
            scaled_toNat_Interactions = [scaled_NatNat_Interactions scaled_NfNat_Interactions];
            scaled_toNf_Interactions = [scaled_NatNf_Interactions scaled_NfNf_Interactions];
           
            % Store the interaction vector
            interactionVectors{row, col} = scaled_toNat_Interactions(:);
            interactionVectors{row, col + n2} = scaled_toNf_Interactions(:);
    
            % Assign growth rate to the growth rate matrix
            Nat_growthRateMatrix(row, col) = muVec(Nat_species);
            Nf_growthRateMatrix(row, col) = muVec(Nf_species);
    
        end
    
    end

end