function species_colormap = makeSpeciesColorMap(Ntot)         

    cmap = jet(Ntot);
    c_jet = colormap(cmap);
    idx = randperm(size(c_jet, 1));
    c_new = c_jet(idx(1:Ntot),:);
    species_colormap = colormap(c_new);

end