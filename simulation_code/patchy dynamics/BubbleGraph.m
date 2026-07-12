function fig = BubbleGraph(n1, n2, V, speciesAssignments, site_colormap, species_colormap, sz, lw)
                                   
    for ix=1:length(V) % ix is site index
                                                                
         site_index = ix;
                                 
        % Calculate row and column indices from linear index
                                    
        row = mod(site_index - 1, n1) + 1;
                                    
        col = ceil(site_index / n1);
                                          
        % Calculate the rotated positions
                                    
        rotated_row = col;  % Swap row and column positions
                                    
        rotated_col = n1 - row + 1;  % Adjust column position
                                
        sID = speciesAssignments(row, col);

        site_color = site_colormap(site_index, :);
        species_color = species_colormap(sID, :);
                              
        scatter(rotated_row, rotated_col, sz * V(ix, 1), site_color, 'filled', 'MarkerEdgeColor', species_color, 'LineWidth', lw);
        hold on;
                                                            
    end
                           
    % Remove all non-essential elements from the figure                     
    axis off;
                        
    %   aspect_ratio = n2 / n1;  % Calculate the aspect ratio
    axis equal;  % Set equal aspect ratio
    %   xlim([0, n1 + 0]);  % Adjust x-axis limits
    %   ylim([0, n2 + 0]);  % Adjust y-axis limits
                        
    hold off

end