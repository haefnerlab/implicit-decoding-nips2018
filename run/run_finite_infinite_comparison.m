%Code to compare finite infinite samples used for inference in orer to
%decode probability of orientation of a grating image

% close all
% clear
N = 21; %Number of neurons
stim_num = N; %Number of stimulus = number of neurons
sig_eb = 0.1; %In NIPS paper \sigma_{exp-brain}
pixel_noise_std = 0.5; %Pixel Noise of images (gratings)
prior = 0.1; %Sparse prior of spiking of neurons
ang_val1 = 0;
ang_val2 = pi;
[G,pix] = tools.PFgenerator(N,ang_val1,ang_val2*(N-1)/N); %generates PFs
params = tools.ModelParams(G,N,pixel_noise_std,prior,pix); %assigns parameters for inference
angle = linspace(ang_val1,ang_val2*(N-1)/N,N);
Im = tools.StimGenerator(params.n_neurons,ang_val1,ang_val2*(params.n_neurons-1)/params.n_neurons);
fix_item = fix(params.n_neurons/2)+1; %Image whose orientation is to be determined
Im_selected = Im(:,fix_item);
[p_inf,mu_vb] = results.InfiniteSampleImplicitCoding_overOrientations(Im,Im_selected,params,sig_eb);


n_samples = [2, 20, 200, 2000];
repeats = 100;
p_fin = zeros(N,length(n_samples));
spk_counts = zeros(N,length(n_samples));
thin_gap = 10;

%% Run the simulation
for i=1:length(n_samples)
    disp(i);
    for j=1:repeats
        Im_selected_noisy = Im(:,fix_item) + sig_eb * randn(params.pixels,1);
        n_total_samples = thin_gap * n_samples(i);
        idx = (1:thin_gap:n_total_samples);
        samples = tools.InferenceGibbs(params, Im_selected_noisy, n_total_samples, 1);
        x_samples = squeeze(samples(:,idx,1));
        p_fin(:,i) = p_fin(:,i) + results.FiniteSampleImplicitCoding_overOrientations(Im,params,sig_eb,x_samples);
        spk_counts(:,i) = spk_counts(:,i) + squeeze(sum(x_samples,2));
    end
    p_fin(:,i) = p_fin(:,i)/repeats;
    spk_counts(:,i) = spk_counts(:,i)/repeats;
end
%% Plot the figure for comparison of probability over decoded abstract variable orientation
color = ['g','y','m','b','c'];
figure(1);
for i=1:length(n_samples)
    p1 = plot(angle,p_fin(:,i),color(i));
    p1.LineWidth = 1;
    hold on;
end
p = plot(angle,p_inf,'k--');
p.LineWidth = 2;
legend({'2 sample', '20 sample', '200 samples', '2000 samples', 'Inf Samples'});
xlabel('Orientations of gratings shown as images')
ylabel('Probabilities over abstract variable (Orientation), given orientated gratings')
title('Finite and infinite case converge for large samples')
ind_mx = find(p_inf==max(max(p_inf),max(p_fin)));
% make sure this is the fix item in the infinite_Sampling code as well
p2 = line([angle(fix_item) angle(fix_item)], [max(p_inf) 0]);
p2.Color = [1 0 0];
p2.LineWidth = 2;
set(gca, 'XTick', sort([angle(fix_item), get(gca, 'XTick')]));
axis('tight')

%% Plot the figure for comparison of spikecounts over decoded abstract variable orientation
figure(2);
for i=1:length(n_samples)
    p1 = plot(angle,spk_counts(:,i),color(i));
    p1.LineWidth = 1;
    hold on;
end
legend({'2 sample', '20 sample', '200 samples', '2000 samples'});
xlabel('Orientations of gratings shown as images')
ylabel('SpikeCounts over Orientation, given orientated gratings')
title('Finite and infinite case converge for large samples')
axis('tight')




