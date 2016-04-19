clear all; clc; close all;
%%
size = 256;
x = linspace(-1,1,size);
[X,Y] = meshgrid(x,-x);
rad = 1;
circ = X.^2+Y.^2 <= rad^2;
rho = sqrt(X.^2+Y.^2);
theta = atan2(Y,X);

field = zeros(size,size);
load('PS7_Zernike_coeffs.mat');

Z = zeros(256,256,25);
Z(:,:,1) = rho.^0;
Z(:,:,2) =   rho.*cos(theta);
Z(:,:,3) =   rho.*sin(theta);
Z(:,:,4) =   2*rho.^2-1;
Z(:,:,5) =   rho.^2.*cos(2*theta);
Z(:,:,6) =   rho.^2.*sin(2*theta);
Z(:,:,7) =   (3*rho.^2-2).*rho.*cos(theta);
Z(:,:,8) =   (3*rho.^2-2).*rho.*sin(theta);
Z(:,:,9) =    6*rho.^4-6*rho.*2+1;
Z(:,:,10) =   rho.^3*cos(3*theta);
Z(:,:,11) =   rho.^3*sin(3*theta);
Z(:,:,12) =   (4*rho.^2-3).*rho.^2.*cos(2*theta);
Z(:,:,13) =   (4*rho.^2-3).*rho.^2.*sin(2*theta);
Z(:,:,14) =   (10*rho.^4-12*rho.^2+3).*rho.*cos(theta);
Z(:,:,15) =   (10*rho.^4-12*rho.^2+3).*rho.*sin(theta);
Z(:,:,16) =   20*rho.^6-30*rho.^4+12*rho.^2-1;
Z(:,:,17) =   rho.^4*cos(4*theta);
Z(:,:,18) =   rho.^4*sin(4*theta);
Z(:,:,19) =   (5*rho.^2-4).*rho.^3.*cos(3*theta);
Z(:,:,20) =  (5*rho.^2-4).*rho.^3.*sin(3*theta);
Z(:,:,21) =   (15*rho.^4-20*rho.^2+6).*rho.^2.*cos(2*theta);
Z(:,:,22) =   (15*rho.^4-20*rho.^2+6).*rho.^2.*sin(2*theta);
Z(:,:,23) =   (35*rho.^6-60*rho.^4+30*rho.^2-4).*rho.*cos(theta);
Z(:,:,24) =   (35*rho.^6-60*rho.^4+30*rho.^2-4).*rho.*sin(theta);
Z(:,:,25) =   70*rho.^8-140*rho.^6+90*rho.^4-20*rho.^2+1;

wvf = 0;
for i = 1:25
    wvf = wvf+Zernike_Coefficients(i)*Z(:,:,i);
end
imshow(wvf.*circ,[]);
axis on
xlabel('\rho_x');
ylabel('\rho_y');
ax = gca; % current axes
ax.FontSize = 12;
% ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
ax.XLim = [-1 1];
ax.YLim = [-1 1];
