function site_colormap = makeSiteColorMap(n1, n2)         

    Z = n1 * n2;

    cmap = jet(Z);
    c_jet = colormap(cmap);
    idx = randperm(size(c_jet, 1));
    c_new = c_jet(idx(1:Z),:);
    site_colormap = colormap(c_new);

end