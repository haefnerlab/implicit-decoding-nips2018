function [G, pix] = PFgenerator(N, start_angle, stop_angle)
% TOOLS.PFGENERATOR generates Gabor shaped projective fields with orientations spanning between the
% start_angle and the stop_angle. 
%
% - N: the number of neurons, i.e number of projective fields to create
% - start_angle, stop_angle: specifies the range to span the orientation of the Gabor like projective fields

angle = linspace(start_angle, stop_angle, N);
X1 = (-5: 0.2: 5);
Y1 = (-5: 0.2: 5);
w = length(X1);
[X, Y] = meshgrid(X1, Y1);
pix = w * w;
G = zeros(pix, N);
for i = 1:N
    Z = tools.rf(2, 1, 1/.56, 0.0, X, Y, angle(i));
    G(:, i) = Z(:);
    G(:, i) = G(:, i) / sqrt(G(:, i)'*G(:, i));
end
end