function [Im] = StimGenerator(StimN, start_ang, stop_ang)
% STIMGENERATOR generate grating stimulus to be used as input to the network of neurons while doing
% inference in the binary sparse coding model with Gabor shaped projective fields of various
% orientations.
%
% Note that we deliberately use even-symmetric gratings with zero phase offset both here and in the
% definition of the model's projective fields so that results will be periodic in pi rather than 2pi
%
% - StimN: number of stimulus to generate
% - start_ang, stop_ang: determines the range of angles to span

angle = linspace(start_ang, stop_ang, StimN);
X1 = (-5:0.2: 5);
Y1 = (-5:0.2:5);
w = length(X1);
[X, Y] = meshgrid(X1, Y1);
pix = w * w;

Im = zeros(pix, StimN);
for i = 1:StimN
    Z = stm(1/.56, X, Y, angle(i));
    Im(:, i) = Z(:);
    Im(:, i) = Im(:, i)/sqrt((Im(:, i)'*Im(:, i)));
end
end

function temp = stm(k, x, y, rotate)
y_t = -x * sin(rotate) + y * cos(rotate);
temp = cos(k*y_t);
end