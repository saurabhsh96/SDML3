%PhiIntegration for uniform current
function Iphi = PhiInt(k0, er, h, ksw, phi, JFT)
    dph = phi(2,1) - phi(1,1);
    Iphi = (sum(((abs(JFT).^2).*(cos(phi)).^2), 'all')).*dph;
end