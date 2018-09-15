%%Code to reproduce the effect of offset as nuisance variable on decoded
%orientation

% close all
% clear all
N = 105;
stim_num = N;
sig_eb = 0.1;
pixel_noise_std = 0.5;
prior = 0.01;
[G,pix] = tools.PFgenerator(N,0,2*pi*(N-1)/N);
params = tools.ModelParams(G,N,pixel_noise_std,prior,pix);
lo = -2;
up = 2;
num = 5;
dc = linspace(lo,up,num);
p = zeros(length(dc),N);
parfor i=1:length(dc)
    disp(i);
    p(i,:) = results.InfiniteSampleImplicitCoding_overOrientations_offset_test(params,sig_eb,dc(i),lo,up);
end

%% Generate figure
offset = linspace(lo,up,num);
angle = linspace(0,2*pi*(N-1)/N,params.n_neurons);
figure();
for i=1:num
    plot(angle,p(i,:),'o-','LineWidth',2)
    hold on;
end
axis('tight')
legend(arrayfun(@(c) ['offset = ' num2str(c)], offset, 'UniformOutput', false));
% make sure this is the fix item in the infinite_Sampling code as well
fix_item = fix(params.n_neurons/2) + 1;
p2 = line([angle(fix_item) angle(fix_item)], [max(p(end,:)) min(p(end,:))]);
p2.Color = [1 0 0]; 
p2.LineWidth = 1;
yticks([]);
set(gca, 'XTick', sort([angle(fix_item), get(gca, 'XTick')]));
