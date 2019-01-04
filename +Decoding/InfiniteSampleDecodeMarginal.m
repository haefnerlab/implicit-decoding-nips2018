function p = InfiniteSampleDecodeMarginal(Templates, params, x_bar, sig_eb)
reconstruction = params.pf * x_bar(:);
dists = sum((Templates - reconstruction).^2) / (sig_eb^2);
dists = dists - min(dists(:));
p = exp(-1/2 * dists);
p = p / sum(p(:));
end