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

Z = zeros(size,size,25);
Z(:,:,1) =   rho.^0;
Z(:,:,2) =   rho.*cos(theta);
Z(:,:,3) =   rho.*sin(theta);
Z(:,:,4) =   2*rho.^2-1;
Z(:,:,5) =   rho.^2.*cos(2*theta);
Z(:,:,6) =   rho.^2.*sin(2*theta);
Z(:,:,7) =   (3*rho.^2-2).*rho.*cos(theta);
Z(:,:,8) =   (3*rho.^2-2).*rho.*sin(theta);
Z(:,:,9) =    6*rho.^4-6*rho.^2+1;
Z(:,:,10) =   rho.^3.*cos(3*theta);
Z(:,:,11) =   rho.^3.*sin(3*theta);
Z(:,:,12) =   (4*rho.^2-3).*rho.^2.*cos(2*theta);
Z(:,:,13) =   (4*rho.^2-3).*rho.^2.*sin(2*theta);
Z(:,:,14) =   (10*rho.^4-12*rho.^2+3).*rho.*cos(theta);
Z(:,:,15) =   (10*rho.^4-12*rho.^2+3).*rho.*sin(theta);
Z(:,:,16) =   20*rho.^6-30*rho.^4+12*rho.^2-1;
Z(:,:,17) =   rho.^4.*cos(4*theta);
Z(:,:,18) =   rho.^4.*sin(4*theta);
Z(:,:,19) =   (5*rho.^2-4).*rho.^3.*cos(3*theta);
Z(:,:,20) =   (5*rho.^2-4).*rho.^3.*sin(3*theta);
Z(:,:,21) =   (15*rho.^4-20*rho.^2+6).*rho.^2.*cos(2*theta);
Z(:,:,22) =   (15*rho.^4-20*rho.^2+6).*rho.^2.*sin(2*theta);
Z(:,:,23) =   (35*rho.^6-60*rho.^4+30*rho.^2-4).*rho.*cos(theta);
Z(:,:,24) =   (35*rho.^6-60*rho.^4+30*rho.^2-4).*rho.*sin(theta);
Z(:,:,25) =   70*rho.^8-140*rho.^6+90*rho.^4-20*rho.^2+1;

wvfsys = 0;
for i = 1:25
    wvfsys = wvfsys+Zernike_Coefficients(i)*Z(:,:,i);
end
wvfsys = wvfsys.*circ;

imshow(wvfsys, [],'XData', [-1 1], 'YData', [-1 1]), colorbar;
colormap(jet)
axis on
xlabel('\rho_x');
ylabel('\rho_y');

dz = -1000:500:1000;
NA = 1/sqrt(1+100);

for i = 1:length(dz)
    wvf = wvfsys+0.5*NA^2.*dz(i).*(X.^2+Y.^2);
    Eobs = fftshift(ifftn(ifftshift(circ.*exp(-1i*2*pi*wvf))));
    figure;
    psf = abs(Eobs).^2;
    imshow(circ.*psf.^0.25,[],'XData', [-1 1], 'YData', [-1 1]), colorbar;
    colormap(jet)
    axis on
    xlabel('\rho_x');
    ylabel('\rho_y');
end

%% c)
for i = 1:length(dz)
    knifeedge = Y >= 0;
    wvf = wvfsys+0.5*NA^2.*dz(i).*(X.^2+Y.^2);
    Eobs = fftshift(fftn(ifftshift(knifeedge.*fftshift(ifftn(ifftshift(circ.*exp(-1i*2*pi*wvf)))))));
    figure;
    psf = abs(Eobs).^2;
    imshow(circ.*psf.^0.25,[],'XData', [-1 1], 'YData', [-1 1]), colorbar;
    colormap(jet)
    axis on
    xlabel('\rho_x');
    ylabel('\rho_y');
end

%% d)
%calculate impulse response
center = fix(size/2)+1;

wvf = wvfsys;
A(:,:,i) = fftshift(ifftn(ifftshift(circ.*exp(-1i*2*pi*wvf))));

target = A(:,:,i);
h = fftshift(ifftn(ifftshift(target)));
h2 = abs(h).^2;
temp = fftshift(fftn(ifftshift(h2)));
mtf = abs(temp/temp(center,center));
imshow(mtf.^0.25, 'xdata', [-1 1], 'ydata', [-1 1]),colorbar;
colormap(jet);
axis on
xlabel('\rho_x [cm]');
ylabel('\rho_y [cm]');

firstzero = min(find(mtf(center, center:size)<=1e-3));
plotx = linspace(0,2,firstzero);
ploty = mtf(center,center:center+firstzero-1);
figure;
plot(plotx, ploty);
xlabel('f_x/f_0');
ylabel('MTF_{fx}');

firstzero = min(find(mtf(center:size, center)<=1e-3));
plotx = linspace(0,2,firstzero);
ploty = mtf(center:center+firstzero-1,center);
figure;
plot(plotx, ploty);
xlabel('f_y/f_0');
ylabel('MTF_{fy}');




