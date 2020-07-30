%Finding Iz
function Iz = Zint(k0, er, h, ksw, mode)
    %Impedences
    zeta0 = 120*pi;
    
    %In Slab
    ks = sqrt(er).*k0;
    zetaS = zeta0./sqrt(er);
    
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
    
    constant = (Zup.*Zdown)./dKgDel;   
    
    TermA = (zetaS./ks).*(abs(1j.*constant./(Zs.*sin(kzs.*h)))).^2.*(h/2)*(1 + sinc((2.*kzs.*h)/pi));
    TermB = (zeta0./k0).*(abs(constant./Z0)).^2.*(1./(2.*sqrt((ksw.^2) - (k0.^2))));
    
    Iz = TermA + TermB;
    %Const term in residue, not dependant on integration
%     const_termRes = ((abs(Zup.*Zdown./dKgDel)).^2).*(zeta0./k0);
%     
%     %Integration 
%     tillH = (h./2).*(1+sinc(2.*kzs.*h./pi))./(((abs(Zs.*sin(kzs.*h))).^2).*er); 
%     aboveH = 1./(((abs(Z0)).^2).*2.*(sqrt(ksw.^2 - k0.^2)));
%     Iz = (tillH + aboveH).*const_termRes;
end