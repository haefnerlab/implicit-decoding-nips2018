% MODELPARAMS initializes the parameters of the binary sparse coding model
% that will be used to perform inference.
% G: the matrix of projective fields or basis functions
% Neurons: number of neurons in the population
% pixel_noise_std: standard deviation of external pixel noise
% prior: sparse prior of spiking of each neuron in the network
% pixels: size of image patches and projective fields/basis of neurons

function params = ModelParams(G, Neurons, pixel_noise_std, prior, pixels)
params = struct('pf', G, 'n_neurons', Neurons, 'pixel_std', pixel_noise_std, ...
    'prior', prior, 'pixels', pixels);
end