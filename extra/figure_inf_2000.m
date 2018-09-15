figure();
load('mu_vb.mat')
load('post_2000.mat')
load('post_vb.mat')
load('samples_2000.mat')
N = 21;
angle = linspace(0,pi*(N-1)/N,N);
fix_item = fix(N/2)+1;
% subplot(2,2,[1,3])
p1 = plot(angle,exp(p_final-logsumexp(p_final)),'-o');
p1.LineWidth = 2;
hold on;
p = plot(angle,p_inf,'k--');
p.LineWidth = 1;
legend({'2000 samples', 'Inf Samples'});
ind_mx = find(p_inf==max(p_inf));
p2 = line([angle(fix_item) angle(fix_item)], [max(p_inf) 0]);
p2.Color = [1 0 0];
p2.LineWidth = 1;
% set(gca, 'XTick', sort([angle(fix_item), get(gca, 'XTick')]));
% axis('tight')
%% FR
figure();
% subplot(2,2,2)
p1 = plot(angle,sum(x_samples,2),'-*');
p1.LineWidth = 2;
hold on;
% subplot(1,2,4)
p2 = plot(angle,mu_vb*2000,'k-*');
p2.LineWidth = 2;
legend({'2000 samples', 'Variational Bayes'});