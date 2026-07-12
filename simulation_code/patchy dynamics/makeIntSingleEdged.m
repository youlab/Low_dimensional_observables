
% % Function to process the existing interactions to 'single edged'
function [a, q] = makeIntSingleEdged(a, q)
         
         [u, v] = find(a > 0);
        
            for i=1:length(u)
        
                a(v(i), u(i)) = 0;
        
            end
        
            [u, v] = find(q > 0);
        
            for i=1:length(u)
        
                q(v(i), u(i)) = 0;
        
            end

end