function p_inf_cont = InfiniteSampleImplicitCoding_overOrientations_contrast_test(params, sig_eb, true_contrast, T)
%INFINITESAMPLEIMPLICITCODING_OVERORIENTATIONS_CONTRAST_TEST returns the decoded probability of an
%abstract variable (here orientation) for a given level of nuisance contrast

Templates = tools.StimGenerator(T, 0, pi*(T-1)/T);
mu_vb = zeros(1,params.n_neurons);
fix_item = round(T/2); % present the 'middle' template to the model
variational_trials = 1000;

% Get _expected_ mean response after 1000 noise draws added to the template at contrast 'con'
for i=1:variational_trials
    Im_noisy = true_contrast * Templates(:,fix_item) + sig_eb * randn(params.pixels,1);
    mu_vb = mu_vb + tools.variational_bayes(params, Im_noisy);
end
mu_vb = mu_vb / variational_trials;
reconstruction = params.pf*mu_vb(:);

% Assume that the decoder does not have access to the true contrast and that contrast itself has a
% uniform prior in [0 1], and marginalize contrast out
n_contrast = 500;
contrast = linspace(0, 1, n_contrast);
prior_contrast = ones(size(contrast)) / n_contrast;
log_prob_c_s = zeros(T, n_contrast);
for j=1:n_contrast
    template_dist = contrast(j) * Templates - reconstruction;
    log_prob_c_s(:,j) = -sum(template_dist.^2 / sig_eb^2) / 2;
end
prob_c_s = exp(log_prob_c_s - max(log_prob_c_s(:)));
p_inf_cont = prob_c_s * prior_contrast(:);
p_inf_cont = p_inf_cont / sum(p_inf_cont);
end

