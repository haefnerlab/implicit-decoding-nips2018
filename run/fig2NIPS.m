%Code to get the decoded proability of orientation for inference done with
%finite samples

% close all
% clear all
N = 21; %Number of neurons
stim_num = N; %Number of stimulus = number of neurons
sig_eb = 0.1; %In NIPS paper \sigma_{exp-brain}
pixel_noise_std = 0.5; %Pixel Noise of images (gratings)
prior = 0.1; %Sparse prior of spiking of neurons
ang_val1 = 0;
ang_val2 = pi;
[G,pix] = tools.PFgenerator(N,ang_val1,ang_val2*(N-1)/N); %generates PFs
params = tools.ModelParams(G,N,pixel_noise_std,prior,pix); %assigns parameters for inference
repeats = 3;
Im = tools.StimGenerator(params.n_neurons,ang_val1,ang_val2*(params.n_neurons-1)/params.n_neurons);
fix_item = fix(params.n_neurons/2)+1; %Image whose orientation is to be determined
Im_selected = Im(:,fix_item);
n_samples_set = [2,20,200,2000];
%Drawing finite independent samples
figure();
thin_gap = 50;
angle = linspace(ang_val1,ang_val2*(N-1)/N,N);
for i=1:length(n_samples_set)
    n_samples = n_samples_set(i);% number of effective samples needed
    n_total_samples = thin_gap*n_samples;
    if (n_samples == max(n_samples_set))
        repeats = 1;
    end
    for j=1:repeats 
        samples = tools.InferenceGibbs(params, Im_selected, n_total_samples, 1);
        idx = (1:thin_gap:n_total_samples);
        idx_selected = datasample(idx,n_samples,'Replace',false);
        x_samples = squeeze(samples(:,idx_selected,1));
        
        
        p_fin = results.FiniteSampleImplicitCoding_overOrientations(Im,params,sig_eb,x_samples);
        
        %% plotting the results
        subplot(2,length(n_samples_set),i)
        if (n_samples == max(n_samples_set))
            [p_inf,mu_vb] = results.InfiniteSampleImplicitCoding_overOrientations(Im,Im_selected,params,sig_eb);
        end
        p = plot(angle,p_fin,'-o');
        p.LineWidth = 2;
        xlabel('Orientation');
        ylabel('Probability');
        axis('tight');
        ylim([0 0.5])
        legend(num2str(n_samples));
        hold on;
        if (n_samples == max(n_samples_set))
            p=plot(angle,p_inf,'k--');
            p.LineWidth = 4;
            legend(num2str(n_samples),'Inf samples');
            hold on;
        end
        
        subplot(2,length(n_samples_set),i+length(n_samples_set))
        s = sum(x_samples,2);
        p = plot(angle,s,'-o');
        p.LineWidth = 2;
        xlabel('Orientation');
        ylabel('Spike counts');
        axis('tight');
        legend(num2str(n_samples));
        hold on;
        if (n_samples == max(n_samples_set))
            p=plot(angle,mu_vb*n_samples,'k--');
            p.LineWidth = 4;
            legend(num2str(n_samples),'VB');
            hold on;
        end
    end
end