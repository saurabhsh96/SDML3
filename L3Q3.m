%Lecture 3, Question 3
clear;
close all;

%% Defining Inputs

%Dimensions
h = 2e-3;

%EM
er = 10;
freq = 1e9:0.1e9:20e9;
c = 3e8;
%lam = c./freq;
%k0 = 2.*pi./lam;

%Impedances
zeta0 = 120*pi;
zetaS = zeta0./(sqrt(er));

%Current orientation
J = [1, 0, 0];

%% Q3-1 Power in FF and SW by ED normalized to the same in Free space

%For space
%Meshgrid in Space
drad = pi/180;
ph = (0:1:360)*pi/180;
th = linspace(eps,89,181)*pi/180;
% th = linspace(eps, pi/2-drad, 90);
% ph = linspace(eps, 2*pi-drad, 360);
[thi, phi] = meshgrid(th, ph);
dth = thi(1, 2) - thi(1, 1);
dph = phi(2, 1) - phi(1, 1);

%For Slab
%rho_vec = 0:0.000001:0.000001; %make sure it is in FF otherwise can create issue?
rho_vec = 1.00;
%rho_vec = 1:0.001:1.1;
%phi_vec = eps:drad/4:2*pi;
phi_vec = (0:0.1:360)*pi/180;
[rho, phiR] = meshgrid(rho_vec, phi_vec);
% dphiR = phiR(2, 1) - phiR(1, 1);
%Dir = zeros([size(thi) size(freq, 2)]);

%Spectral %FS %SW %Elem
PradElem = zeros(size(freq));
PradFElem = zeros(size(freq));
PradSWElem = zeros(size(freq));

%Spectral %FS %SW %Elem
PradUnif = zeros(size(freq));
%PradFUnif = zeros(size(freq));
PradSWUnif = zeros(size(freq));

%Spectral %FS %SW %Elem
PradPWS = zeros(size(freq));
%PradFPWS = zeros(size(freq));
PradSWPWS = zeros(size(freq));

%Efficiencies
% etaElem = zeros(size(freq));
% etaUnif = zeros(size(freq));
% etaPWS = zeros(size(freq));

%Ksw taken from L3Q1 Loading only TM0 data
kSW = load('kSW.mat');

%% Loop
for ind = 1:size(freq, 2)
    %% Common variables
    fr = freq(ind);
    lam = c/fr;
    k0 = 2*pi/lam;
    ks = k0*sqrt(er);
    ksw = kSW.kgRTM(ind);
    
    keq = (k0+ks)./2;
    
%     if(isnan(ksw))
%         continue;
%     end
    
    kxs = k0.*sin(thi).*cos(phi);
    kys = k0.*sin(thi).*sin(phi);
    kzs = k0.*cos(thi);
    kRho = sqrt(kxs.^2 + kys.^2); 

    %Observation point
    r = 1;
    z = r.*cos(thi);
    z_dash = h;
    %z = r;
    %z_dash = 0;

    %Tx-line equivalent circuits
    [vTM, vTE, iTM, iTE] = trxline_GroundSlab(k0, er, h, zeta0, zetaS, kRho, z);

    %Green's function calculations
    [Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SGFejF(k0, er, ...
        kxs, kys, vTM, vTE, iTM, iTE, zeta0, zetaS);
    
    %% Prad in FS by spectral and only element in FS by elementary
    %For elementary current
    JFT = ones([size(J, 2) size(kxs)]);
    JFT(2, :, :) = 0;
    JFT(3, :, :) = 0;
    
    %JFT = 1;
    %Calculating Electric field
    [~, ~, Emag, ~] = Field(k0, kzs, r, thi, phi, Gxx, Gyx,...
        Gzx, Gxy, Gyy, Gzy, JFT, z, z_dash);

    %Spectral version of directivity function Dir(:,:,ind)
    [~, PradElem(1, ind)] = DirectivityF(Emag, 1, r, thi, dth, dph);
    
    %Free space -> Assuming no stratified medium at all, not even the
    %ground plane, As taking till theta = 90 deg, prad is doubled in below
    %calculations DirF(:,:,ind)
    [~, PradFElem(1, ind)] = Directivity(fr, 0, 0, 1, r, thi, phi);
    
    JFT = ones([size(J, 2) size(phiR)]);
    JFT(2, :, :) = 0;
    JFT(3, :, :) = 0;
    %Calculation of SW power
    PradSWElem(1, ind) = PswTM(k0, er, h, ksw, 0, 0, squeeze(JFT(1,:,:)), phiR, "TM");
    
    %% Prad calculations for PWS
    %JFT of the current
    L = 5.3e-3;
    W = 0.5e-3;
    
    %For PWS
    JFT = CurrentFT(keq, kxs, kys, L, W, J);
    
    %Calculating Electric field
    [Eth, Eph, Emag, ~] = Field(k0, kzs, r, thi, phi, Gxx, Gyx,...
        Gzx, Gxy, Gyy, Gzy, JFT, z, z_dash);

    %Spectral version of directivity function Dir(:,:,ind)
    [~, PradPWS(1, ind)] = DirectivityF(Emag, 1, r, thi, dth, dph);
    
    %Free space -> Assuming no stratified medium at all, not even the
    %ground plane, As taking till theta = 90 deg, prad is doubled in below
    %calculations DirF(:,:,ind)
    %[~, PradFPWS(1, ind)] = Directivity(fr, L, W, 1, r, thi, phi);
    
    %For SW
    kx = ksw.*cos(phiR);
    ky = ksw.*sin(phiR);
    
    %keq = k0.*sqrt((1+er)./2);
    
    JFT = CurrentFT(keq, kx, ky, L, W, J);
 
    %SW power
    PradSWPWS(1, ind) = (PswTM(k0, er, h, ksw, L, W, squeeze(JFT(1,:,:))', phiR, "TM"));
    %% Calculation for unifrom source
    %Psw = (1/2)(ksw^2/2pi)IzIphi
    %JFT of the current
    L = 25e-3;
    W = 25e-3;
    
    %keq = (k0+ks)./2;
    
    %For PWS
    JFT = CurrentFTUnif(keq, kxs, kys, L, W, J);
    
    %Calculating Electric field
    [~, ~, Emag, ~] = Field(k0, kzs, r, thi, phi, Gxx, Gyx,...
        Gzx, Gxy, Gyy, Gzy, JFT, z, z_dash);

    %Spectral version of directivity function Dir(:,:,ind)
    [~, PradUnif(1, ind)] = DirectivityF(Emag, 1, r, thi, dth, dph);
    
    %Free space -> Assuming no stratified medium at all, not even the
    %ground plane, As taking till theta = 90 deg, prad is doubled in below
    %calculations DirF(:,:,ind)
    %[~, PradFUnif(1, ind)] = Directivity(fr, L, W, 1, r, thi, phi);
    
    kx = ksw.*cos(phiR);
    ky = ksw.*sin(phiR);
    
    %keq = k0.*sqrt((1+er)./2);
    
    JFT = CurrentFTUnif(keq, kx, ky, L, W, J);
    
    %SW power
    PradSWUnif(1, ind) = (PswTM(k0, er, h, ksw, L, W, squeeze(JFT(1,:,:))', phiR, "TM"));
end

%% Plotting
figure(); %Q3-1
PradReq = PradElem./(2*PradFElem); %FS
PradReq1 = PradSWElem(1,2:length(PradSWElem))./(2*PradFElem(1,2:length(PradSWElem)));
%Temporarily removing NaN -> Think of something later
plot(freq./(10^9), pow2db(PradReq), 'LineWidth', 1.5, 'DisplayName', "Prad"); hold on;
plot(freq(1,2:length(PradSWElem))./(10^9), pow2db(PradReq1), 'LineWidth', 1.5, 'DisplayName', "Psw");
title('P_{rad} (Normalized w.r.t. P_{rad} in FS) vs. Frequency');
xlabel('Frequency (in GHz)');
ylabel('Normalized P_{rad} (in dB) ');
legend show;
grid on;

%% Efficiency calculations
etaElem = PradElem./(PradElem + PradSWElem);
etaUnif = PradUnif./(PradUnif + PradSWUnif);
etaPWS = PradPWS./(PradPWS + PradSWPWS);

%Plotting efficiencies
figure(); %Q3-1
plot(freq./(10^9), etaElem, 'LineWidth', 1.5, 'DisplayName', "Elementary"); hold on;
plot(freq./(10^9), etaUnif, 'LineWidth', 1.5, 'DisplayName', "Uniform");
plot(freq./(10^9), etaPWS, 'LineWidth', 1.5, 'DisplayName', "PWS");
%ylim([0.6, 1]);
title('Efficiency vs. Frequency');
xlabel('Frequency (in GHz)');
ylabel('Efficiency');
legend show;
grid on;
