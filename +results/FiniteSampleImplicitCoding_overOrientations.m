
%FINITESAMPLINGIMPLICITCODING_OVERORIENTATIONS is a function that returns the
%decoded probability of an abstract variable of interest, in this case
%orientation of an image shown based on sampled responses of a network of
%neurons in a binary sparse coding model.
%params: paramteres to simulate the results in the NIPS paper
%sig_eb: in the paper describes as \sigma_{exp-brain}; i.e sigma noise that
%emerges when moving from the experimenter to the brain interface
%x_samples: the number of samples drawn in inference over binary sparse
%coding model
%Returns probability of orientation of the image(grating) shown p_inf, for
%infinite number of samples drawn to do inference

function p_fin = FiniteSampleImplicitCoding_overOrientations(Templates,params,sig_eb,x_samples)
size_x_samples = size(x_samples);
n_samples = size_x_samples(2);
var_final = ((sig_eb^2) + params.pixel_std^2/n_samples);
p = zeros(params.n_neurons,1);
for i=1:size(Templates, 2)
    p(i) = sum(log(normpdf(Templates(:,i),params.pf*mean(squeeze(x_samples),2),sqrt(var_final))));
end
p_fin = exp(p - tools.logsumexp(p));
end


