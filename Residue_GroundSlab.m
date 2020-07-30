%Residue function for ground slab
function [VtmR, ItmR] = Residue_GroundSlab(k0, er, h, ksw, z, mode)  
    %Zeta
    zeta0 = 120*pi;
    zetaS = zeta0./sqrt(er);
   
    %In Slab
    ks = sqrt(er).*k0;
    
    %kZ
    kz0 = -1j*sqrt(-((k0^2)-(ksw.^2)));
    kzs = -1j*sqrt(-((ks^2)-(ksw.^2)));

    %Defining according to the mode
    if(mode == "TE")
        Z0 = (zeta0.*k0)./kz0;
        Zs = (zetaS.*ks)./kzs;
        Zup = Z0;
        Zdown = 1j*Zs.*tan(kzs.*h);
    else 
        Z0 = (zeta0.*kz0)./k0;
        Zs = (zetaS.*kzs)./ks;
        Zup = Z0;
        Zdown = 1j*Zs.*tan(kzs.*h);
    end 
    
    %delK
    dK = k0./500;
    
    %Derivative
    dKgDelP = Den_GroundSlab(k0, er, h, (ksw + dK/2), zeta0, mode);
    dKgDelM = Den_GroundSlab(k0, er, h, (ksw - dK/2), zeta0, mode);
    dKgDel = (dKgDelP - dKgDelM)./dK;
   
    %Finding residues
    if(z<=h)
        VtmR = Zup.*Zdown.*sin(kzs.*z)./(dKgDel.*sin(kzs*h));
        ItmR = Zup.*Zdown.*1j.*cos(kzs.*z)./(dKgDel.*sin(kzs*h).*Zs);
    else
        VtmR = Zup.*Zdown.*exp(1j.*kz0.*h).*exp(-1j.*kz0.*z)./(dKgDel);
        ItmR = Zup.*Zdown.*exp(1j.*kz0.*h).*exp(-1j.*kz0.*z)./(dKgDel.*Z0);
    end 
end