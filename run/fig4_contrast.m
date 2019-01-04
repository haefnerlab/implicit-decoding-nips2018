%% Code to reproduce the effect of contrast as nuisance variable on decoded orientation

N = 72;          % Number of neurons
T = 201;         % Number of templates along the stimulus manifold
stim_num = N;    % Number of stimuli = number of neurons
sig_eb = 0.1;    % Pixel noise added to templates. In NeurIPS paper, \sigma_{exp-brain}
sig_model = 0.5; % Pixel noise assumed by the model (we allow for some model-mismatch here)
prior = 0.1;     % Sparse prior of spiking of neurons

[G, pix] = tools.PFgenerator(N, 0, pi*(N-1)/N); % generates PFs
G = G - mean(G, 1); % Ensure that there is no 'DC' component in the projective fields
params = BinarySparseCoding.ModelParams(G, N, sig_model, prior, pix);


contrast = linspace(0, 1, 5);
p_inf_cont = zeros(length(contrast), T);
for i=1:length(contrast)
    fprintf('Running contrast level %d/%d\n', i, length(contrast));
    p_inf_cont(i, :) = Decoding.ContrastNuisance(params, sig_eb, contrast(i), T);
end

%% Generating the figure
angle = linspace(0, pi, T+1);
figure();
colors = repmat(linspace(.8, 0, length(contrast))', 1, 3);
for i = 1:length(contrast)
    plot(angle, [p_inf_cont(i, :) p_inf_cont(i, 1)], '-', 'LineWidth', 2, 'Color', colors(i, :))
    hold on;
end
axis('tight')
legend(arrayfun(@(c) ['contrast = ' num2str(c)], contrast, 'UniformOutput', false));
% make sure this is the fix item in the infinite_Sampling code as well
yticks([]);
set(gca, 'fontsize', 10, 'fontweight', 'bold', 'xtick', [0 pi/4 pi/2 3*pi/4 pi], 'xticklabel', {0, '\pi/4', '\pi/2', '3\pi/4', '\pi'})
