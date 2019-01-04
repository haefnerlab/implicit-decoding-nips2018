function gabor = rf(sig_x, sig_y, k, phi, x, y, rotate)
x_t = x * cos(rotate) + y * sin(rotate);
y_t = -x * sin(rotate) + y * cos(rotate);
magnitude = 1.0/(2*pi*sig_x*sig_y);
window = exp((-((x_t.*x_t)/(2*sig_x^2))-((y_t.*y_t)/(2*sig_y^2))));
gabor = magnitude * window .* cos(k*y_t-phi);
end