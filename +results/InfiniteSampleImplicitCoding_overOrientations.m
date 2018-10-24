function [p_inf,mu_vb] = InfiniteSampleImplicitCoding_overOrientations(Templates,Im_noisy,params,sig_eb)
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
mu_vb = tools.variational_bayes(params,Im_noisy);
p = zeros(params.n_neurons,1);
for i=1:size(Templates, 2)
    p(i) = sum(log(normpdf(Templates(:,i),params.pf*mu_vb(:),(sig_eb))));   
end
p_inf = exp(p - tools.logsumexp(p));
end


