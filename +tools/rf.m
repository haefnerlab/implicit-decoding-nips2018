%RF is a sub function used to generate Gabor shaped projective fields of
%different orientations
function temp =  rf(sig_x,sig_y,k,phi,x,y,rotate)
        temp=(1.0/(2*pi*sig_x*sig_y));
        x_t = x * cos(rotate) + y * sin(rotate);
        y_t = -x * sin(rotate) + y * cos(rotate);
        temp = temp * exp((-((x_t.*x_t)/(2*sig_x^2))-((y_t.*y_t)/(2*sig_y^2))));
        temp = temp.*cos(k*y_t-phi);
        
end