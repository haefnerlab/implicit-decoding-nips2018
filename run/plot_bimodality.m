%Code to make bimodal image of the NIPS paper

% close all
% clear
N = 101; %Number of neurons
stim_num = N; %Number of stimulus = number of neurons
sig_eb = 0.05; %In NIPS paper \sigma_{exp-brain}
pixel_noise_std = 0.5; %Pixel Noise of images (gratings)
prior = 0.01; %Sparse prior of spiking of neurons
[G,pix] = tools.PFgenerator(N,0,pi*(N-1)/N); %generates PFs
params = tools.ModelParams(G,N,pixel_noise_std,prior,pix); %assigns parameters for inference
Im = tools.StimGenerator(params.n_neurons,0.0,pi*(params.n_neurons-1)/params.n_neurons);
lag = 29;
fix_item1 = fix(params.n_neurons/2)+lag;
fix_item2 = fix(params.n_neurons/2)-lag;
Im_selected = Im(:,fix_item1)+Im(:,fix_item2);
Im_selected = Im_selected/sqrt(Im_selected'*Im_selected);
p_inf = zeros(params.n_neurons,1);
repeats = 1;
for i=1:repeats
    %disp(i);
    p_inf = p_inf + results.InfiniteSampleImplicitCoding_overOrientations(Im,Im_selected,params,sig_eb);
end
p_inf = p_inf/repeats;
%% generating figure
angle = linspace(0,2*pi*(N-1)/N,N);
p = plot(angle,p_inf,'o-');
p.LineWidth = 4;
xlabel('Orientations of gratings shown as images')
ylabel('Probabilities over abstract variable (Orientation), given orientated gratings')
title('Finite and infinite case converge for large samples')
ind_mx = find(p_inf==max(p_inf));
% make sure this is the fix item in the infinite_Sampling code as well
p2 = line([angle(fix_item1) angle(fix_item1)], [max(p_inf) min(p_inf)]);
p2.Color = [1 0 0]; 
p2.LineWidth = 2;
set(gca, 'XTick', sort([angle(fix_item1), get(gca, 'XTick')]));
hold on;
p2 = line([angle(fix_item2) angle(fix_item2)], [max(p_inf) min(p_inf)]);
p2.Color = [1 0 0]; 
p2.LineWidth = 2;
set(gca, 'XTick', sort([angle(fix_item2), get(gca, 'XTick')]));
axis('tight')