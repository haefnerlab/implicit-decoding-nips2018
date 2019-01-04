N = 36;          % Number of neurons
T = 201;         % Number of templates along the stimulus manifold
n_repeats = 100; % Number of runs to average over
stim_num = N;    % Number of stimuli = number of neurons
sig_eb = 0.1;    % Pixel noise added to templates. In NeurIPS paper, \sigma_{exp-brain}
sig_model = 0.5; % Pixel noise assumed by the model (we allow for some model-mismatch here)
prior = 0.1;     % Sparse prior of spiking of neurons
dbg = true;      % Flag for debugging mode

% Create binary sparse coding model
[G, pix] = tools.PFgenerator(N, 0, pi*(N-1)/N);
params = BinarySparseCoding.ModelParams(G, N, sig_model, prior, pix);

% Create image templates (gratings)
Templates = tools.StimGenerator(T, 0, pi*(T-1)/T);

% Choose the middle template as the true image
fix_item = round(T/2);
true_template = Templates(:, fix_item);

% Create finite samples plots for each of the following numbers of samples
n_samples_set = [1 10 100];

% Gibbs parameters – ensure the chain is sufficiently burnt in, and thin the samples so they are
% independent (as required by the derivations in the paper)
thin_gap = 10;
burn_in = 500;

n_samples = n_samples_set(end);                    % max number of effective samples needed
n_total_samples = burn_in + thin_gap * n_samples;  % total number of samples to draw

%% Draw samples

if dbg, disp('Drawing samples for finite-sample cases'); end

samples = zeros(N, n_total_samples, n_repeats);

for j = 1:n_repeats
    fprintf('Repeat %02d/%02d\n', j, n_repeats);

    % Add noise to the template image with variance set by sig_eb
    noisy_image = true_template + sig_eb * randn(params.pixels, 1);
    
    % Run Gibbs sampling
    samples(:, :, j) = BinarySparseCoding.InferenceGibbs(params, noisy_image, n_total_samples, 1);
end

thin_samples = samples(:, burn_in+1:thin_gap:n_total_samples, :);

%% If debugging, check autocorrelation

if dbg
    dbg_fig = figure;
    hold on;
    for j=1:10
        acf = [];
        for n=1:N
            [this_acf, lags] = autocorr(squeeze(samples(n, burn_in+1:end, j)), 'NumLags', 5*thin_gap);
            acf = vertcat(acf, this_acf);
        end
        mean_acf = nanmean(acf, 1);
        stderr_acf = nanstd(acf, [], 1) / sqrt(N);
        h(j) = errorbar(lags, mean_acf, stderr_acf);
    end
    yl = ylim;
    plot([thin_gap thin_gap], yl, '--r');
    text(thin_gap + 1, 0.1, sprintf('thinning = %d', thin_gap));
    title('Gibbs sampler autocorrelation');
    legend(h, arrayfun(@(j) sprintf('run %02d', j), 1:10, 'UniformOutput', false));
    pause;
    try
        close(dbg_fig);
    catch
    end
end

%% VB solution as an approximation for infinite-time limit

if dbg, disp('Computing VB solution'); end

[p_inf, mu_vb] = Decoding.InfiniteSampleVBApproximation(Templates, true_template, params, sig_eb);

%% Infinite (really, very large finite) solution

if dbg, disp('Drawing 10k samples for approximately-infinite solution'); end

inf_sample_approx = 10000;
n_total_samples = burn_in + inf_sample_approx;
inf_samples = BinarySparseCoding.InferenceGibbs(params, true_template, n_total_samples, 1);
x_bar = squeeze(mean(inf_samples(:, burn_in+1:end, :), 2));
p_fin = Decoding.InfiniteSampleDecodeMarginal(Templates, params, x_bar, sig_eb);

%% Before plotting, select which 3 runs to plot to give best visualizations

% Find runs where the first sample had only a single active neuron (easier to visualize)
run_ids = find(sum(thin_samples(:, 1, :), 1) == 1);
cell_ids = arrayfun(@(i) find(thin_samples(:, 1, i)), run_ids);
[~, cellsort] = sort(cell_ids);
run_ids = run_ids(cellsort);
% Out of those runs, after sorting by the spiking neuron's preferred orientation, select those at
% the 25th, 50th, and 75th percentile to get a good visualizable spread of initial values.
ids = run_ids([round(length(run_ids)/4) round(length(run_ids)/2) round(3*length(run_ids)/4)]);

%% Compute finite-sample decoded distributions

if dbg, disp('Decoding finite-sample distributions'); end

% NOTE: we assume here that the integral in equation (5) of the paper is constant. This is
% reasonable in this case where we've used entirely rotationally-symmetric stimuli and projective
% fields, but will not be true in general. This means that the 'finite sample' decoder is identical
% to the infinite sample decoder.

ps = zeros(n_repeats, T, length(n_samples_set));

for j=1:n_repeats
    for i=1:length(n_samples_set)
        x_bar = mean(thin_samples(:, 1:n_samples_set(i), j), 2);
        ps(j, :, i) = Decoding.InfiniteSampleDecodeMarginal(Templates, params, x_bar, sig_eb);
    end
end

%% Plot everything

if dbg, disp('Plotting results'); end

axformat = {'fontsize', 10, 'fontweight', 'bold', 'xtick', [0 pi/4 pi/2 3*pi/4 pi], 'xticklabel', {0, '\pi/4', '\pi/2', '3\pi/4', '\pi'}};
cols = hsv(3);

% Create x-axis indices
angle_T = linspace(0, pi*(T-1)/T, T);
angle_N = linspace(0, pi*(N-1)/N, N);

% Plot finite-samples, individual runs
for j=1:3
    for i=1:length(n_samples_set)
        subplot(2, length(n_samples_set)+1, i)
        hold on
        plot(angle_T, ps(ids(j), :, i), '-', 'linewidth', 2, 'Color', cols(j, :));
        set(gca, axformat{:})
        xlabel('Stim. Ori. (s)');
        ylim([0 0.15]);
        
        subplot(2, length(n_samples_set)+1, 1+i+length(n_samples_set))
        hold on
        stem(angle_N, squeeze(sum(thin_samples(:, 1:n_samples_set(i), ids(j)), 2)), '.', 'Color', cols(j, :), 'linewidth', 1);
        xlabel('Pref. Ori.');
        set(gca, axformat{:})
    end
end

% Plot finite-samples average run
for i=1:length(n_samples_set)
    subplot(2, length(n_samples_set)+1, i)
    hold on
    plot(angle_T, mean(ps(:, :, i), 1), 'k-', 'linewidth', 2);
    
    subplot(2, length(n_samples_set)+1, 1+i+length(n_samples_set))
    hold on
    plot(angle_N, squeeze(mean(sum(thin_samples(:, 1:n_samples_set(i), :), 2), 3)), 'k-', 'linewidth', 2);
end

% Plot infinite-samples case
subplot(2, length(n_samples_set)+1, length(n_samples_set)+1)
hold on
plot(angle_T, p_inf, 'k', 'linewidth', 2);
plot(angle_T, p_fin, 'linewidth', 2, 'Color', [.5 .5 .5]);
plot(angle_T, p_inf, 'k', 'linewidth', 2);
xlabel('Stim. Ori. (s)');
set(gca, axformat{:})
ylim([0 0.08]);

subplot(2, length(n_samples_set)+1, 2*length(n_samples_set)+2)
hold on
stem(angle_N, sum(inf_samples, 2), '.', 'Color', [.5 .5 .5], 'linewidth', 1);
plot(angle_N, mu_vb*inf_sample_approx, 'k-', 'linewidth', 2);
xlabel('Pref. Ori.');
set(gca, axformat{:})

subplot(2, length(n_samples_set)+1, 1)
ylabel('p(s|spikes)');
subplot(2, length(n_samples_set)+1, length(n_samples_set)+2)
ylabel('Spike Counts');
set(gcf, 'color', 'white')