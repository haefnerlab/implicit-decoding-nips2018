%INFINITESAMPLEIMPLICITCODING_OVERORIENTATIONS_OFFSET_TEST is a function
%that returns the decoded probability of an abstract variable here
%orientation for different levels of added offset of the grating image
%shown.

function p_inf_offset = InfiniteSampleImplicitCoding_overOrientations_offset_test(params,sig_eb,dc_fix,lo,up)
Im = tools.StimGenerator(params.n_neurons,0.0,2*pi*(params.n_neurons-1)/params.n_neurons);
mu_vb = zeros(1,params.n_neurons);
% angle = linspace(0,pi*(params.n_neurons-1)/params.n_neurons,params.n_neurons);
fix_item = fix(params.n_neurons/2)+1;
variational_trials = 1000;
for i=1:variational_trials
    Im_noisy = Im(:,fix_item) + dc_fix +  sig_eb * randn(params.pixels,1);
    mu_vb = mu_vb + tools.variational_bayes(params,Im_noisy);
end
mu_vb = mu_vb/variational_trials;
p = zeros(params.n_neurons,1);
p_inf_offset = zeros(params.n_neurons,1);
dc = unifrnd(lo,up,1,100);
for j=1:100
    for i=1:params.n_neurons
        Im1 = Im(:,i)+ones(params.pixels,1)*dc(j);
        p(i) = sum(log(normpdf(Im1,params.pf*mu_vb(:),(sig_eb))));  
    end
    p_inf_offset = p_inf_offset + exp(p - tools.logsumexp(p));
end
p_inf_offset = p_inf_offset/sum(p_inf_offset);
end


