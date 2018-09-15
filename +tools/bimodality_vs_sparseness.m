%% Test bimodality as a function of sparseness

xs = linspace(-3, 3, 101);
[xx, yy] = meshgrid(xs);
betas = [2 1 .5]; % try priors of the form exp(|x|^beta) for different betas (2=gaussian, 1=exponential)
mu_x = [.5 .5];
cov_x = [.5 -.3; -.3 .5];

dxy = [xx(:), yy(:)] - mu_x;
dxy_scaled = dxy / cov_x;
log_likelihood = reshape(-sum(dxy_scaled .* dxy, 2), size(xx));
likelihood = log_to_pdf(log_likelihood);

figure;
contour(xx, yy, likelihood, linspace(0, max(likelihood(:)), 7));
title('likelihood');
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
axis square;

figure;
for ib = 1:length(betas)
    log_priors{ib} = -abs(xx).^betas(ib) -abs(yy).^betas(ib);
    prior_pdf = log_to_pdf(log_priors{ib});
    
    log_posteriors{ib} = log_priors{ib} + log_likelihood;
    post_pdf = log_to_pdf(log_posteriors{ib});
    
    subplot(2, length(betas), ib);
    hold on;
    contour(xx, yy, prior_pdf, linspace(0, max(prior_pdf(:)), 7));
    title(['Prior (\beta = ' num2str(betas(ib)) ')']);
    set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
    axis square;
    
    subplot(2, length(betas), length(betas) + ib);
    hold on;
    contour(xx, yy, post_pdf, linspace(0, max(post_pdf(:)), 7));
    title(['Posterior (\beta = ' num2str(betas(ib)) ')']);
    set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
    axis square;
end

%% Helper
function p = log_to_pdf(log_p)
p = exp(log_p - max(log_p(:)));
end