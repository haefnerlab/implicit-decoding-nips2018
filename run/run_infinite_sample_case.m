%Code to get the decoded proability of orientation for inference done with
%infinite samples

% close all
% clear
N = 101; %Number of neurons
stim_num = N; %Number of stimulus = number of neurons
sig_eb = 0.1; %In NIPS paper \sigma_{exp-brain}
pixel_noise_std = 0.5; %Pixel Noise of images (gratings)
prior = 0.1; %Sparse prior of spiking of neurons
[G,pix] = tools.PFgenerator(N,-ang_val,ang_val*(N-1)/N); %generates PFs
ang_val = pi/2;
params = tools.ModelParams(G,N,pixel_noise_std,prior,pix); %assigns parameters for inference

Im = tools.StimGenerator(params.n_neurons,-ang_val,ang_val*(params.n_neurons-1)/params.n_neurons);
fix_item = fix(params.n_neurons/2)+1; %Image whose orientation is to be determined
Im_selected = Im(:,fix_item);
[p_inf,mu_vb] = results.InfiniteSampleImplicitCoding_overOrientations(Im,Im_selected,params,sig_eb);

%% generating figure
angle = linspace(-ang_val,ang_val*(N-1)/N,N);
p = plot(angle,p_inf);
p.LineWidth = 2;
xlabel('Neuron Orientation');
ylabel('Probability of orientation of selected image Im selected');
axis('tight')