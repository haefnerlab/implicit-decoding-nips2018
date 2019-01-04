function [p_inf, mu_vb] = InfiniteSampleVBApproximation(Templates, Im_noisy, params, sig_eb)
% DECODING.INFINITESAMPLEVBAPPROXIMATION returns the decoded probability of an abstract variable of
% interest that defines the manifold of Templates. Uses the variational bayes approximation to
% marginal probabilities in the infinite-samples limit.
%
% Templates: [pixels x T] densely-sampled templates evenly spaced along the manifold of interest
% Im_noisy: Image on which inference is performed, with noise already added
% params: see BinarySparseCoding.ModelParams
% sig_eb: true noise in the image --> brain process (\sigma_{exp --> brain} in the paper)
%
% Returns p_inf, a [1 x T] probability mass function over the provided templates, as well as mu_vb,
% the mean-field VB solution

mu_vb = BinarySparseCoding.VariationalBayes(params, Im_noisy);
p_inf = Decoding.InfiniteSampleDecodeMarginal(Templates, params, mu_vb, sig_eb);
end
