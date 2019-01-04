function x_samples = InferenceGibbs(params, image, n_samples, n_repeats)
% BINARYSPARSECODING.INFERENCEGIBBS do inference in the binary sparse coding model using gibbs
% sampling. Samples are initialized from the prior. Returns (neurons x n_samples x n_repeats) array
% of sampled values (0s and 1s).
%
% x_samples = BINARYSPARSECODING.INFERENCEGIBBS(params, image, samples, repeats)

% Flatten the image
image = image(:);

feedforward = image' * params.pf;
R = -params.pf' * params.pf;
pix_var = params.pixel_std^2;
log_prior = log(params.prior) - log(1 - params.prior);

x_samples = zeros(params.n_neurons, n_samples, n_repeats);

assert(numel(image) == params.pixels);

% No outer loop over repeats (vectorized)

% Sample 1st value from the prior
x_samples(:, 1, :) = rand(params.n_neurons, n_repeats) < params.prior;

% Loop sequentially over samples
for s = 2:n_samples
    % Copy state s-1 to state s
    x_samples(:, s, :) = x_samples(:, s-1, :);
    
    % Loop random permutation over cells (note: this is the same random permutation per repeat due
    % to vectorization, but different for each sample)
    for k = randperm(params.n_neurons)
        % Sample x_k
        idx_other = [1:k-1 k+1:params.n_neurons];
        drive = (feedforward(k) + R(k, k) / 2 + squeeze(x_samples(idx_other, s, :))' * R(idx_other, k)) / pix_var;
        prob = tools.sigmoid_x(log_prior + drive);
        x_samples(k, s, :) = rand(size(prob)) < prob;
    end
end

end