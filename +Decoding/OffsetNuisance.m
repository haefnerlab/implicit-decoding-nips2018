function p_inf_offset = OffsetNuisance(params, sig_eb, true_off, lo, up, T)
% INFINITESAMPLEIMPLICITCODING_OVERORIENTATIONS_OFFSET_TEST returns the decoded probability of an
% abstract variable (here orientation) for a given level of added luminance offset

Templates = tools.StimGenerator(T, 0, pi*(T-1)/T);
Templates = Templates - mean(Templates, 1); % Ensure that there is no 'DC' component in the templates themselves

mu_vb = zeros(1, params.n_neurons);
fix_item = round(T/2); % present the 'middle' template to the model
variational_trials = 1000;

% Get _expected_ mean response (approximated with variational bayes) after 1000 noise draws added to
% the template at offset 'true_off'
for i = 1:variational_trials
    Im_noisy = true_off + Templates(:, fix_item) + sig_eb * randn(params.pixels, 1);
    mu_vb = mu_vb + BinarySparseCoding.VariationalBayes(params, Im_noisy);
end
mu_vb = mu_vb / variational_trials;
reconstruction = params.pf*mu_vb(:);

% Assume that the decoder does not have access to the true offset and that the offset itself has a
% uniform prior in [lo up], and marginalize it out
n_offset = 500;
offsets = linspace(lo, up, n_offset);
prior_contrast = ones(size(offsets)) / n_offset;
log_prob_off_s = zeros(T, n_offset);
for j = 1:n_offset
    template_dist = offsets(j) + Templates - reconstruction;
    log_prob_off_s(:, j) = -sum(template_dist.^2 / sig_eb^2) / 2;
end
prob_off_s = exp(log_prob_off_s - max(log_prob_off_s(:)));
p_inf_offset = prob_off_s * prior_contrast(:);
p_inf_offset = p_inf_offset / sum(p_inf_offset);
end


