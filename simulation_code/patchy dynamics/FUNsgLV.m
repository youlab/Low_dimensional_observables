function dy = FUNsgLV(t, y, delta, F, G, GR)
    
    dy =  GR .* y .* (1 - y) - (delta .* y ./ (F*y + 1) ) - (G*y .* y) ;

end