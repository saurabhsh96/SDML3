%SW fields
function [Erho, Ez, Hphi] = TMSwFields(k0, ksw, er, VtmR, ItmR, rho, phi, z, J, L, W, h)
    %zeta
    zeta0 = 120*pi;
    zetaS = zeta0./sqrt(er);
    ks = k0.*sqrt(er);
    
    %Constant terms in EFF
    C = 1j.*sqrt(ksw./(2*pi)).*exp(1j.*pi./4);
    const_term = C.*cos(phi).*exp(-1j.*ksw.*rho)./sqrt(rho);
    
    %JFT
    %kx = ksw.*cos(phi);
    %ky = ksw.*sin(phi);
    
    %keq = k0.*sqrt((1+er)./2);
    %JFT = CurrentFT(keq, kx, ky, L, W, J);
    
    %Jx = squeeze(JFT(1,:,:));
    Jx = 1;
    
    %EFF calculations
    Erho = VtmR.*Jx.*const_term;
    Hphi = ItmR.*Jx.*const_term;
    
    if(z<=h)
        Ez = -((zetaS.*ksw)./ks).*ItmR.*Jx.*const_term;
    else 
        Ez = -((zeta0.*ksw)./k0).*ItmR.*Jx.*const_term;
    end 
end