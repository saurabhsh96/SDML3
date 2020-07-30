%lecture 3, Q1
clear;
close all;

%% Input Definitions

%Material
er = [5,10];
h = 2e-3;

%EM
freq = 1e9:0.1e9:20e9;
c = 3e8;
lam = c./freq;
k0 = 2.*pi./lam;

%Wave impedances
eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;
%zeta0 = (sqrt(mu_0/(eps_0*1)));
zeta0 = 120*pi;

%Meshgrid
% drad = pi/180;
% th = linspace(eps, pi/2-drad, 90);
% ph = linspace(eps, 2*pi-drad, 360);
% [thi, phi] = meshgrid(th, ph);
% dth = thi(1, 2) - thi(1, 1);
% dph = phi(2, 1) - phi(1, 1);

%% Finding required

%kRho for TM
freqL = length(freq);
%Er in question
erQ = er(1);
ks = k0(freqL)*sqrt(erQ);
kRho = k0(freqL):0.1:ks; %for larger

%D value for largest freq
DlargestTM = Den_GroundSlab(k0(freqL), erQ, h, kRho, zeta0, "TM");
DlargestTE = Den_GroundSlab(k0(freqL), erQ, h, kRho, zeta0, "TE");

%Required plots
figure();
plot(kRho./k0(freqL), abs(DlargestTM), 'LineWidth', 1.5, 'DisplayName', "TM"); 
hold on;
plot(kRho./k0(freqL), abs(DlargestTE), 'LineWidth', 1.5, 'DisplayName', "TE");
title(['Denominator function to obtain guess poles (\epsilon_r = ', num2str(erQ), ')']);
%title('Initial Guess, for highest freq');
xlabel('k^g_\rho/k_0');
ylabel('D(k^g_\rho)');
grid on;
legend show;
%ylim([0, 10]);

%Finding initial value
[~, Itm] = min(abs(DlargestTM));
[~, Ite] = min(abs(DlargestTE));

%Required KRho value
kgItm = kRho(Itm);
kgIte = kRho(Ite);

%NR method
kgPrevTM = kgItm;
kgPrevTE = kgIte;

%Defining required kgRho
kgRTM = zeros(size(freq));
kgRTE = zeros(size(freq));

%Intital defination
kgRTM(freqL) = kgItm;
kgRTE(freqL) = kgIte;

%Freq loop
for ind = size(freq, 2):-1:2
    
    %Normalized kg
    kgnTM = kgPrevTM./k0(ind);
    kgnTE = kgPrevTE./k0(ind);
    
    %Kg for next step
    kgTM = kgnTM.*k0(ind-1);
    kgTE = kgnTE.*k0(ind-1);
    
    %Required kg
    kgRTM(ind - 1) = findprop(k0(ind-1), erQ, h, kgTM, zeta0, "TM");
    kgRTE(ind - 1) = findprop(k0(ind-1), erQ, h, kgTE, zeta0, "TE");
    
    %Next guess
    kgPrevTM = kgRTM(ind - 1);
    kgPrevTE = kgRTE(ind - 1);
    
    %Only evaluate real values
    if(imag(kgRTM(ind - 1))~=0)
        kgRTM(ind - 1) = 0;
    end
    if(imag(kgRTE(ind - 1))~=0)
        kgRTE(ind - 1) = 0;
    end
end

%Plotting Krho vs. freq
figure();
plot(freq./10^9, real(kgRTM./k0), 'LineWidth', 1.5, 'DisplayName', "TM"); 
hold on;
plot(freq./10^9, real(kgRTE./k0), 'LineWidth', 1.5, 'DisplayName', "TE");
title(['Guess pole point vs. Freq at \epsilon_r = ', num2str(erQ)]);
ylabel('k^g_\rho/k_0');
xlabel('Frequency (in GHz)');
ylim([1, 3]);
grid on;
legend show;