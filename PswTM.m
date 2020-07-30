%PradSWTM elementary current
function Psw = PswTM(k0, er, h, ksw, L, W, JFT, phi, mode)
    const_term = (1./2).*((ksw.^2)./(2.*pi));
    Iz = Zint(k0, er, h, ksw, mode);
    
    if(L == 0 && W == 0)
        Iphi = pi;
    else 
        Iphi = PhiInt(k0, er, h, ksw, phi, JFT);
    end 
    
    Psw = const_term.*Iz.*Iphi;
end