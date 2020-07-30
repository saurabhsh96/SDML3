%Lecture 3, Q2: TM EFF, surface waves
clear;
close all;

%% Defining Input

%Dimension
h = 2e-3;

%EM
er = 10;
freq = 10e9;
c = 3e8;
lam = c/freq;
k0 = 2*pi/lam;

%Unnecesary for elementary source
L = lam./200000;
W = lam./200000;

%% Q2-1 and 2-3: 
%Real and imaginary part variation of EF as a function of radial distance
%and w.r.t phi

%Current
J =[1, 0, 0];

%Krho taken from previous question TM0 and freq = 10, er = 10, in Q1, gives
%the required krho
kRho = 264.7765;

%Defining coordinates
rho_vec = 1:0.001:1.1; %make sure it is in FF otherwise can create issue?
drad = pi/180;
phi_vec = eps:drad:2*pi;
[rho, phi] = meshgrid(rho_vec, phi_vec);

%Resudue calculations
ksw = kRho;
z = 1e-3;
[VtmR, ItmR] = Residue_GroundSlab(k0, er, h, ksw, z, "TM");

%Field calculations %Careful about L and W
[Erho, Ez, Hphi] = TMSwFields(k0, ksw, er, VtmR, ItmR, rho, phi, z, J, L, W, h);

%Magnitudes
EMag = sqrt((abs(Erho)).^2 + (abs(Ez)).^2);
Emax = max(max(EMag(1,:,:)));

%Plotting the Erho and Ez when phi = 0; z = h/2
figure();
plot(rho(1,:), real(Erho(1,:)), 'LineWidth', 1.2, 'DisplayName', 'Re(E_\rho)'); hold on;
plot(rho(1,:), imag(Erho(1,:)), '--', 'LineWidth', 1.2, 'DisplayName', 'Im(E_\rho)');
plot(rho(1,:), real(Ez(1,:)), 'LineWidth', 1.2, 'DisplayName', 'Re(E_z)');
plot(rho(1,:), imag(Ez(1,:)), '--', 'LineWidth', 1.2, 'DisplayName', 'Im(E_z)');
xlabel('\rho [in m]');
ylabel('E_\rho, E_z');
title('Variation of E_\rho, E_z w.r.t. radial distance');
%xlim([1, 2]);
legend show;
grid on;
hold off;

%Plotting for constant rho and z and varying phi 
figure();
plot(phi(:,1)./drad, abs(Erho(:,1)./Emax), 'LineWidth', 1.2, 'DisplayName', '|(E_\rho)|'); hold on;
plot(phi(:,1)./drad, abs(Ez(:,1)./Emax), 'LineWidth', 1.2, 'DisplayName', '|(E_z)|');
plot(phi(:,1)./drad, abs(EMag(:,1)./Emax), 'LineWidth', 1.2, 'DisplayName', '|(E_{tot})|');
title('Amplitude variation of the electric field w.r.t. \phi');
ylabel('Normalized E_\rho and E_z');
xlabel('\phi [in deg]');
legend show;
grid on;
hold off;

%Polar plots
figure();
polarplot(phi(:,1), abs(Erho(:,1)./Emax), 'LineWidth', 1.2, 'DisplayName', '|(E_\rho)|'); hold on;
polarplot(phi(:,1), abs(Ez(:,1)./Emax), 'LineWidth', 1.2, 'DisplayName', '|(E_z)|');
polarplot(phi(:,1), abs(EMag(:,1)./Emax), 'LineWidth', 1.2, 'DisplayName', '|(E_{tot})|');
title('Polar plot of Amplitude variation electric field w.r.t. \phi');
legend show;
hold off;

%% Q2-2 Plotting against Z

zVar = linspace(eps, 5*h, 100);

%Varying with z
Erho1 = zeros([size(zVar,2) size(rho)]);
Ez1 = zeros([size(zVar,2) size(rho)]);
Hphi1 = zeros([size(zVar,2) size(rho)]);

for ind = 1:size(zVar, 2)
    %Resudue calculations
    [VtmR1, ItmR1] = Residue_GroundSlab(k0, er, h, ksw, zVar(ind), "TM");

    %Field calculations
    [Erho1(ind,:,:), Ez1(ind,:,:), Hphi1(ind,:,:)] = TMSwFields(k0, ksw, ...
        er, VtmR1, ItmR1, rho, phi, zVar(ind), J, L, W, h);
end

%Plotting the Erho and Ez when phi = 0, rho = 1m;
figure();
%Required rho index
%reqRhoInd = find(ismembertol(rho(1,:), 0.1));
%Plotting
EMag = sqrt((abs(Erho1)).^2 + (abs(Ez1)).^2);
Emax1 = max(max(max(EMag)));
plot(zVar.*1e3, abs(squeeze(Erho1(:,1,1))./Emax1), 'LineWidth', 1.5, 'DisplayName', '|E_\rho|'); hold on;
%plot(zVar.*10e3, imag(squeeze(Erho1(:,1,reqRhoInd))), '--', 'LineWidth', 1.2, 'DisplayName', 'Im(E_\rho)');
plot(zVar.*1e3, abs(squeeze(Ez1(:,1,1))./Emax1), 'LineWidth', 1.5, 'DisplayName', '|E_z|');
plot(zVar.*1e3, abs(squeeze(EMag(:,1,1))./Emax1), 'LineWidth', 1.5, 'DisplayName', 'E_{tot}');
xlabel('z [in mm]');
ylabel('E_\rho, E_z Normalized w.r.t E_{tot}');
title('Variation of E_\rho, E_z w.r.t. Z');
%xlim([1, 2]);
legend show;
grid on;
hold off;