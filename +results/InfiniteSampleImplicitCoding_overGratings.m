%INFINITESAMPLINGIMPLICITCODING_OVERGRATINGS is a function that returns the
%decoded probability of an abstract variable of interest, in this case
%orientation of an image shown based on sampled responses of a network of
%neurons in a binary sparse coding model.
%params: paramteres to simulate the results in the NIPS paper
%sig_eb: in the paper describes as \sigma_{exp-brain}; i.e sigma noise that
%emerges when moving from the experimenter to the brain interface
%Returns probability of orientation of the image(grating) shown p_inf, for
%infinite number of samples drawn to do inference

function p_inf = InfiniteSampleImplicitCoding_overGratings(params,sig_eb)
Im = tools.StimGenerator(params.n_neurons);
mu_vb = zeros(1,params.n_neurons);
angle = linspace(0,pi,params.n_neurons);
fix_item = fix(params.n_neurons/2)+1;
variational_trials = 1000;
for i=1:variational_trials
    Im_noisy = Im(:,fix_item) +  sig_eb * randn(params.pixels,1);
    mu_vb = mu_vb + tools.variational_bayes(params,Im_noisy);
end
mu_vb = mu_vb/variational_trials;
p = zeros(params.n_neurons,1);
figure();
for i=1:params.n_neurons
    p(i) = sum(log(normpdf(Im(:,i),params.pf*mu_vb(:),(sig_eb))));

end
p_inf = exp(p - tools.logsumexp(p));
ind_mx = find(p_inf==max(p_inf));
plot(angle,p_inf,'o-')
hold on;
plot((1:1:params.n_neurons),p_inf,'o-')
plot(angle(ind_mx),p_inf(ind_mx),'ro');
xlabel('Orientations of Gabor Projective Fields')
ylabel('Probability of abstract variable (Orientation)')
title(['Infinite samples: plot for grating image of orientation ',num2str(angle(fix_item))]);
axis('tight')
p = line([angle(fix_item) angle(fix_item)], [max(p_inf) 0]);
p.Color = [1 0 0]; 
set(gca, 'XTick', sort([angle(fix_item), get(gca, 'XTick')]));
title(['Infinite samples: plot for grating image of orientation ',num2str(angle(fix_item))]);% plot(angle(ind_mx),p1(ind_mx),'ro');
end
