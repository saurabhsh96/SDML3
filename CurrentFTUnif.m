% Routine to calculate the Fourier Transform of a current
function [currFT] = CurrentFTUnif(k0, kl, kw, L, W, J)
    %kl represents the orientation of Length of dipole
    %kw represents the orientation of Width of dipole
    %Assuming Uniform as current Dist.
    %Assuming I0 = 1
    %J = [1; 0; 0];
    Tx = sinc(kw*W/2/pi);
    Lx = sinc(kl*L/2/pi);
    
    currFT(1,:,:) = L.*Lx.*Tx*J(1);
    currFT(2,:,:) = L.*Lx.*Tx*J(2);
    currFT(3,:,:) = L.*Lx.*Tx*J(3);
end