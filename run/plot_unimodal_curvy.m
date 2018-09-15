%Code to make unimodal image of the NIPS paper

% close all
% clear all
N = 105;
stim_num = N;
sig_eb = 0.1;
pixel_noise_std = 0.5;
prior = 0.01;
[G,pix] = tools.PFgenerator(N,-2*pi,2*pi*(N-1)/N);
params = tools.ModelParams(G,N,pixel_noise_std,prior,pix);
Im = tools.StimGenerator(params.n_neurons,0.0,pi*(params.n_neurons-1)/params.n_neurons);
lag = 35;
fix_item1 = fix(params.n_neurons/2)-lag;
angle = linspace(-2*pi,2*pi*(N-1)/N,N);
indxs = [fix_item1+lag:5:length(angle)];
Im_selected = Im(:,fix_item1);
for i=1:length(indxs)
    Im_selected = Im_selected + 0.25*Im(:,indxs(i));
end
p_inf = zeros(params.n_neurons,1);
repeats = 1;
for i=1:repeats
    %disp(i);
    p_inf = p_inf + results.InfiniteSampleImplicitCoding_overOrientations(Im,Im_selected,params,sig_eb);
end
p_inf = p_inf/repeats;
%% generating figure
p_new = p_inf(10:40);
p_new(end+1:end+31)= p_inf(10:40);
p_new(end+1:end+length(p_inf(10:end)))=p_inf(10:end);
ind_mx = find(p_new==max(p_new));
p1 = plot([1:1:length(p_new)],p_new,'o-','LineWidth',4);
p1.Color = [0 0 1];
% make sure this is the fix item in the infinite_Sampling code as well
p2 = line([ind_mx ind_mx], [max(p_new) min(p_new)]);
p2.Color = [1 0 0]; 
p2.LineWidth = 2;
