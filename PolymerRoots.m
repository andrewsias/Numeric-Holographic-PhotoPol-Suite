%==========================================================================
% HologramHarmonics: Model hologram development via the harmonic expansion
% method.  All equations derived in Mathematica file "Fick's Law solution
% wiht sinusoidal forcing.nb" For FRL.
%
% Aug 2019 RRM
%
% The solution to the coupled first order DE's utilizes the roots of a
% characteristic polynomial.  Solve this only once for efficiency.
% Gm = D_m / (Lambda^2 R_p) where D_m and R_p are the diffusivity and peak
% reaction rate of the monomer;
% Vr = visibility of radical distribution
%==========================================================================
function Proots = PolymerRoots(Gm,Vr)
    
    Poly4 = 1;
    Poly3 = 8+224.*Gm.*pi.^2;
    Poly2 = 24+1344.*Gm.*pi.^2+12544.*Gm.^2.*pi.^4+(-4).*Vr.^2;
    Poly1 = 32+2688.*Gm.*pi.^2+50176.*Gm.^2.*pi.^4+147456.*Gm.^3.*pi.^6+(-16) ...
        .*Vr.^2+(-576).*Gm.*pi.^2.*Vr.^2;
    Poly0 = 2.*(8+147456.*Gm.^3.*pi.^6+(-8).*Vr.^2+Vr.^4+(-64).*Gm.*pi.^2.*(( ...
        -14)+9.*Vr.^2)+(-512).*Gm.^2.*pi.^4.*((-49)+18.*Vr.^2));
    
    Proots = roots([Poly4 Poly3 Poly2 Poly1 Poly0]);
    
end
    
   
      
      