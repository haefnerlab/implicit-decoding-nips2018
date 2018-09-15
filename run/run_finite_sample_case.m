%Code to get the decoded proability of orientation for inference done with
%finite samples

% close all
% clear all
N = 101; %Number of neurons
stim_num = N; %Number of stimulus = number of neurons
sig_eb = 0.1; %In NIPS paper \sigma_{exp-brain}
pixel_noise_std = 0.5; %Pixel Noise of images (gratings)
prior = 0.1; %Sparse prior of spiking of neurons
[G,pix] = tools.PFgenerator(N,-ang_val,ang_val*(N-1)/N); %generates PFs
params = tools.ModelParams(G,N,pixel_noise_std,prior,pix); %assigns parameters for inference
ang_val = pi/2;
Im = tools.StimGenerator(params.n_neurons,-ang_val,ang_val*(params.n_neurons-1)/params.n_neurons);
fix_item = fix(params.n_neurons/2)+1; %Image whose orientation is to be determined
Im_selected = Im(:,fix_item);

%Drawing finite independent samples
thin_gap = 50;
n_samples = 2;% number of effective samples needed
n_total_samples = thin_gap*n_samples;

samples = tools.InferenceGibbs(params, Im_selected, n_total_samples, 1);
idx = (1:thin_gap:n_total_samples);
idx_selected = datasample(idx,n_samples,'Replace',false);
x_samples = squeeze(samples(:,idx_selected,1));

%% Runs the implicit decoding code with finite samples drawn during inference

p_fin = results.FiniteSampleImplicitCoding_overOrientations(Im,params,sig_eb,x_samples);

%% plotting the results
figure();
angle = linspace(-ang_val,ang_val*(N-1)/N,N);
p = plot(angle,p_fin,'-o');
p.LineWidth = 2;
xlabel('Neuron Orientation');
ylabel('Probability of orientation of selected image Im selected');
axis('tight');