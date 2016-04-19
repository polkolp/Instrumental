clear all; clc; close all;
size = 256;
%%
for i = 1:4
    circ = zeros(size,size);
    x = (1:size) - fix(size/2)-1;
    center = find(x==0);
    [X,Y] = meshgrid(x,-x);
    %aperture radius 256/64;
    rad = fix(size/8);
    circ = X.^2+Y.^2 <= rad^2;
    obs = X.^2+Y.^2 <= (i-1)*rad;
    pupil = circ-obs;
    A(:,:,i) = pupil;
end

%calculate impulse response
for i = 1:4
    target = A(:,:,i);
    h = fftshift(ifftn(ifftshift(target)));
    h2 = abs(h).^2;
    temp = fftshift(fftn(ifftshift(h2)));
    mtf = abs(temp/temp(center,center));
    firstzero = min(find(mtf(center, center:size)<=1e-3));
    plotx = linspace(0,2,firstzero);
    ploty = mtf(center,center:center+firstzero-1);
    plot(plotx, ploty);
    hold on
end
xlabel('f_x/f_0');
ylabel('MTF');
legend('obs = 0', 'obs = 0.25\rho', 'obs = 0.5\rho', 'obs = 0.75\rho')
hold off