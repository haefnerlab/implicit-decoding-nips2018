%Code to get the decoded proability of orientation for inference done with
%finite samples

% close all
% clear all
rng('default')
N = 21; %Number of neurons
stim_num = N; %Number of stimulus = number of neurons
sig_eb = 0.1; %In NIPS paper \sigma_{exp-brain}
pixel_noise_std = 0.5; %Pixel Noise of images (gratings)
prior = 0.1; %Sparse prior of spiking of neurons
ang_val1 = 0;
ang_val2 = pi;
[G,pix] = tools.PFgenerator(N,ang_val1,ang_val2*(N-1)/N); %generates PFs
params = tools.ModelParams(G,N,pixel_noise_std,prior,pix); %assigns parameters for inference
repeats = 500;
Im = tools.StimGenerator(params.n_neurons,ang_val1,ang_val2*(params.n_neurons-1)/params.n_neurons);
fix_item = fix(params.n_neurons/2)+1; %Image whose orientation is to be determined
Im_selected = Im(:,fix_item);
n_samples_set = [1,10,100];
figure();
thin_gap = 10;
burn_in=500;
cols=eye(3);
angle = linspace(ang_val1,ang_val2*(N-1)/N,N);

n_samples = n_samples_set(end);% number of effective samples needed
n_total_samples = burn_in+thin_gap*n_samples;

for j=1:repeats
    j
    samples = tools.InferenceGibbs(params, Im_selected, n_total_samples, 1);
    idx = (burn_in+1:thin_gap:n_total_samples);
    idx_selected = datasample(idx,n_samples,'Replace',false);
    x_samples = squeeze(samples(:,idx_selected,1));
    
    
    
    
    %% plotting the results
    
    %subplot(1,2,1)
    %figure(1)
    % if (n_samples == max(n_samples_set))
    [p_inf,mu_vb] = results.InfiniteSampleImplicitCoding_overOrientations(Im,Im_selected,params,sig_eb);
    % end
    for i=1:length(n_samples_set)
        
        
        p_fin = results.FiniteSampleImplicitCoding_overOrientations(Im,params,sig_eb,x_samples(:,1:n_samples_set(i)));
        ps{i}(j,:)=p_fin;
        ss{i}(j,:)=sum(x_samples(:,1:n_samples_set(i)),2);
        
    end
    
    pause(0.1);
end

inf_sample_approx=5000;
n_total_samples = burn_in+thin_gap*inf_sample_approx;
samples = tools.InferenceGibbs(params, Im_selected, n_total_samples, 1);
idx = (burn_in+1:thin_gap:n_total_samples);
idx_selected = datasample(idx,inf_sample_approx,'Replace',false);
x_samples = squeeze(samples(:,idx_selected,1));

p_fin = results.FiniteSampleImplicitCoding_overOrientations(Im,params,sig_eb,x_samples);


alp=0.5;
cols1=(1-alp)*cols+(alp)*ones(3,3);
ids=[9,17,29];
for j=1:3
    for i=1:length(n_samples_set)
        subplot(2,length(n_samples_set)+1,i)
        hold on
        plot(angle,ps{i}(ids(j),:),'o-','linewidth',3,'Color',cols1(j,:));
        set(gca,'fontsize', 20,'fontweight','bold')
        subplot(2,length(n_samples_set)+1,1+i+length(n_samples_set))
        hold on
        plot(angle,ss{i}(ids(j),:),'o-','linewidth',3,'Color',cols1(j,:));
        xlabel('Orientation(s)');
        set(gca,'fontsize', 20,'fontweight','bold')
        
    end
end
for i=1:length(n_samples_set)
    subplot(2,length(n_samples_set)+1,i)
    hold on
    plot(angle,mean(ps{i}),'ko-','linewidth',3);
    subplot(2,length(n_samples_set)+1,1+i+length(n_samples_set))
    hold on
    plot(angle,mean(ss{i}),'ko-','linewidth',3);
    set(gca,'fontsize', 20,'fontweight','bold')
end




subplot(2,length(n_samples_set)+1,length(n_samples_set)+1)
hold on
plot(angle,p_inf,'k--','linewidth',3);
plot(angle,p_fin,'linewidth',3);
set(gca,'fontsize', 20,'fontweight','bold')

subplot(2,length(n_samples_set)+1,2*length(n_samples_set)+2)
hold on
plot(angle,mu_vb*inf_sample_approx,'k--','linewidth',3);
plot(angle,sum(x_samples,2),'linewidth',3);
xlabel('Orientation(s)');
set(gca,'fontsize', 20,'fontweight','bold')


subplot(2,length(n_samples_set)+1,1)
ylabel('Probability');
subplot(2,length(n_samples_set)+1,length(n_samples_set)+2)
ylabel('Spike Counts');
set(gcf,'color','white')