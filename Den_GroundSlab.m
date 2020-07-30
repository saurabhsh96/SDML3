%Ground slab denominator function
function D = Den_GroundSlab(k0, er, h, kRho, zeta0, mode)
    
    %In Slab
    ks = sqrt(er).*k0;
    zetaS = zeta0./sqrt(er);
    
    %kZ
    kz0 = -1j*sqrt(-((k0^2)-(kRho.^2)));
    kzs = -1j*sqrt(-((ks^2)-(kRho.^2)));

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
    
    %Finding Dispersion equation
    D = Zup + Zdown;
end