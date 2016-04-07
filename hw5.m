clear all; clc; close all;
fm = 3e9;
beta = 5;
lambda = 1.053e-6;
c = 3e8;
nu_0 = c/lambda;
%%
peaks = -10:10;
freq = peaks*fm;
smod = besselj(peaks,beta).^2;
bar(freq, smod, 0.1);
xlabel('\nu-\nu_0 [Hz]');
ylabel('S_mod(\nu)/I_0');

%%
samplerate = 256+1;
z = linspace(-0.6035/2,0.6035/2,samplerate);
temp = zeros(length(peaks),length(z));
for i = 1:length(peaks)
    temp(i,:) = besselj(peaks(i),beta).^2*exp(-1i*2*pi*(nu_0-peaks(i)*fm)*2*z/c);
end

gammamod = sum(temp);
vtmod = abs(gammamod);
vtunmod = ones(1,length(z));
figure;
plot(z,vtmod);
hold on
plot(z, vtunmod,'--','Linewidth',1.5);
hold off
xlabel('z [m]');
ylabel('V_T');
legend('modulated','unmodulated');

%%
gammahet = gammamod.*exp(1i*2*pi*nu_0*2*z/c);
recovspec = fftshift(fft(ifftshift(gammahet)));
freqrange = linspace(-6.36e10/2,6.36e10/2,samplerate);
figure();
plot(freqrange,abs(recovspec));
xlabel('\nu-\nu_0 [Hz]');
ylabel('Recovered S_mod(\nu)/I_0');
