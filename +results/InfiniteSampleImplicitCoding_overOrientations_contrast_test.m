%INFINITESAMPLEIMPLICITCODING_OVERORIENTATIONS_CONTRAST_TEST is a function
%that returns the decoded probability of an abstract variable here
%orientation for different levels of contrast of the grating image shown

function p_inf_cont = InfiniteSampleImplicitCoding_overOrientations_contrast_test(params,sig_eb,con)
Im = tools.StimGenerator(params.n_neurons,0.0,2*pi*(params.n_neurons-1)/params.n_neurons);
mu_vb = zeros(1,params.n_neurons);
fix_item = fix(params.n_neurons/2)+1;
variational_trials = 1000;
for i=1:variational_trials
    Im_noisy = con * Im(:,fix_item) +  sig_eb * randn(params.pixels,1);
    mu_vb = mu_vb + tools.variational_bayes(params,Im_noisy);
end
mu_vb = mu_vb/variational_trials;
p = zeros(params.n_neurons,1);
p_inf_cont = zeros(params.n_neurons,1);
cont = unifrnd(0,1,1,100);
for j=1:100
    for i=1:params.n_neurons
      p(i) = sum(log(normpdf(cont(j)* Im(:,i),params.pf*mu_vb(:),(sig_eb))));
    end
    p_inf_cont = p_inf_cont + exp(p - tools.logsumexp(p));
end
p_inf_cont = p_inf_cont/sum(p_inf_cont);
end
