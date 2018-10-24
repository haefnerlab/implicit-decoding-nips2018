%%Code to reproduce the effect of contrast as nuisance variable on decoded
%orientation
% close all
% clear
N = 101; %Number of neurons
stim_num = N; %Number of stimulus = number of neurons
sig_eb = 0.05; %In NIPS paper \sigma_{exp-brain}
pixel_noise_std = 1.0; %Pixel Noise of images (gratings)
prior = 0.01; %Sparse prior of spiking of neurons
[G,pix] = tools.PFgenerator(N,0,pi*(N-1)/N); %generates PFs
params = tools.ModelParams(G,N,pixel_noise_std,prior,pix); %assigns parameters for inference

con = linspace(0,1,5);
p_inf_cont = zeros(length(con),N);
for i=1:length(con)
    disp(i);
    p_inf_cont(i,:) = results.InfiniteSampleImplicitCoding_overOrientations_contrast_test(params,sig_eb,con(i));
end
%% Generating the figure
angle = linspace(0,pi*(N-1)/N,params.n_neurons);
figure();
colors = repmat(linspace(.8, 0, length(con))', 1, 3);
for i=1:length(con)
    plot(angle,p_inf_cont(i,:),'o-','LineWidth',2,'Color',colors(i,:))
    hold on;
end
axis('tight')
legend(arrayfun(@(c) ['contrast = ' num2str(c)], con, 'UniformOutput', false));
% make sure this is the fix item in the infinite_Sampling code as well
fix_item = fix(params.n_neurons/2) + 1;
p2 = line([angle(fix_item) angle(fix_item)], [max(p_inf_cont(end,:)) min(p_inf_cont(end,:))]);
p2.Color = [1 0 0]; 
p2.LineWidth = 1;
yticks([]);
set(gca, 'XTick', sort([angle(fix_item), get(gca, 'XTick')]));
