%% Code to reproduce the effect of contrast as nuisance variable on decoded orientation

N = 32; % Number of neurons
T = 200; % Number of grating orientations (discretization of s)
sig_eb = 0.1; %In NIPS paper \sigma_{exp-brain}
pixel_noise_std = 0.5; %Pixel Noise of images (gratings)
prior = 0.1; %Sparse prior of spiking of neurons
[G, pix] = tools.PFgenerator(N, 0,pi*(N-1)/N); %generates PFs
G = G - mean(G, 1); % Ensure that none of the projective fields have a 'DC' component
params = tools.ModelParams(G, N, pixel_noise_std, prior, pix);

contrast = linspace(0,1,5);
p_inf_cont = zeros(length(contrast),T);
for i=1:length(contrast)
    disp(i);
    p_inf_cont(i,:) = results.InfiniteSampleImplicitCoding_overOrientations_contrast_test(params, sig_eb, contrast(i), T);
end

%% Generating the figure
angle = linspace(0,pi,T+1);
figure();
colors = repmat(linspace(.8, 0, length(contrast))', 1, 3);
for i=1:length(contrast)
    plot(angle,[p_inf_cont(i,:) p_inf_cont(i,1)],'-','LineWidth',2,'Color',colors(i,:))
    hold on;
end
axis('tight')
legend(arrayfun(@(c) ['contrast = ' num2str(c)], contrast, 'UniformOutput', false));
% make sure this is the fix item in the infinite_Sampling code as well
yticks([]);
set(gca, 'fontsize', 10,'fontweight','bold','xtick',[0 pi/4 pi/2 3*pi/4 pi],'xticklabel',{0, '\pi/4', '\pi/2', '3\pi/4', '\pi'})
