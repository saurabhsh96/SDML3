function [Jx] = FTCurrent(k0,er,kx,ky,l,w)

    ks = sqrt(er)*k0;
    keq = (k0+ks)/2;

    T = sinc((ky.*w)./2./pi);
    L = 2*keq.*(cos(kx.*l/2) - cos(keq.*l/2))./((keq.^2 - kx.^2).*sin(keq.*l/2));

    Jx = L.*T;
end
