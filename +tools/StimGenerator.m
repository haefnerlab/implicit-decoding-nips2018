%STIMGENERATOR is a function that generates grating stimulus to be used as
%input to the network of neurons while doing inference in the binary sparse
%coding model with Gabor shaped projective fields of various orientations
%StimN: number of stimulus to generate
%start_ang,stop_ang: determines the range of angles to span

function [Im] = StimGenerator(StimN,start_ang,stop_ang)

% angle array size should be of the size of StimN
angle = linspace(start_ang,stop_ang,StimN);
X1 = (-5:0.2: 5);
Y1 = (-5:0.2:5);
w = length(X1);
[X, Y] = meshgrid(X1, Y1);
pix = w * w;

% create oriented-gabor projective fields
Im = zeros(pix, StimN);
%  figure(200);
for i=1:StimN
    Z = stm(2,1,1/.56,0.0,X,Y,angle(i));
    Im(:,i) = Z(:);
    Im(:,i) = Im(:,i)/sqrt((Im(:,i)'*Im(:,i)));
%     subplot(1, StimN, i);
%     imagesc(reshape(Im(:,i), w, w));
%     colormap('gray');
%     axis 'off';
%     axis image;
end
end

function temp = stm(sig_x,sig_y,k,phi,x,y,rotate)
y_t = -x * sin(rotate) + y * cos(rotate);
temp = sin(k*y_t);
end