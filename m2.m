%==========================================================================
% HologramHarmonics: Model hologram development via the harmonic expansion
% method.  All equations derived in Mathematica file "Fick's Law solution
% wiht sinusoidal forcing.nb" For FRL.
%
% Aug 2019 RRM
%
% Gm = D_m / (Lambda^2 R_p) where D_m and R_p are the diffusivity and peak
% reaction rate of the monomer;
% Vr = visibility of radical distribution = 1/(1 + (2 pi)^2 Gr) where 
% Gr = D_r / (Lambda^2 R_r) where D_r and R_r are the diffusivity and peak
% reaction rate of the radicals;
% tau = (exposure t) R_p
%==========================================================================
function m2val = m2(Gm,Vr,tau,Proots)

m2val = (-Vr.^2/2)*(f(Gm,Vr,tau,Proots(1))+f(Gm,Vr,tau,Proots(2))+...
                    f(Gm,Vr,tau,Proots(3))+f(Gm,Vr,tau,Proots(4)));

end

function fval = f(Gm,Vr,tau,x)

fval = (2.*exp(1).^((1/4).*tau.*x)+144.*exp(1).^((1/4).*tau.*x ...
  ).*Gm.*pi.^2+exp(1).^((1/4).*tau.*x).*x).*((-8)+( ...
  -672).*Gm.*pi.^2+(-12544).*Gm.^2.*pi.^4+(-36864).*Gm.^3.*pi.^6+4.* ...
  Vr.^2+144.*Gm.*pi.^2.*Vr.^2+(-12).*x+(-672).*Gm.*pi.^2.* ...
 x+(-6272).*Gm.^2.*pi.^4.*x+2.*Vr.^2.*x+(-6).* ...
 x.^2+(-168).*Gm.*pi.^2.*x.^2+(-1).*x.^3).^(-1);
  
end