%INFINITESAMPLINGIMPLICITCODING_OVERORIENTATIONS is a function that returns the
%decoded probability of an abstract variable of interest, in this case
%orientation of an image shown based on sampled responses of a network of
%neurons in a binary sparse coding model.
%params: paramteres to simulate the results in the NIPS paper
%sig_eb: in the paper describes as \sigma_{exp-brain}; i.e sigma noise that
%emerges when moving from the experimenter to the brain interface
%Im_selected: Image on which inference is perfomred, i.e whose orientation is to be
%detremined
%Returns probability of orientation of the selected image(grating) shown p_inf, for
%infinite number of samples drawn to do inference and solution of
%variational bayes for the same model mu_vb

function [p_inf,mu_vb] = InfiniteSampleImplicitCoding_overOrientations(Im,Im_selected,params,sig_eb)
mu_vb = zeros(1,params.n_neurons);
variational_trials = 1000;
for i=1:variational_trials
    Im_noisy = Im_selected +  sig_eb * randn(params.pixels,1);
    mu_vb = mu_vb + tools.variational_bayes(params,Im_noisy);
end
mu_vb = mu_vb/variational_trials;
p = zeros(params.n_neurons,1);
for i=1:params.n_neurons
    p(i) = sum(log(normpdf(Im(:,i),params.pf*mu_vb(:),(sig_eb))));   
end
p_inf = exp(p - tools.logsumexp(p));
end


