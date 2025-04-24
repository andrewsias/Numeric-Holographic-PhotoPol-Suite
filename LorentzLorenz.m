% ***********************************************************************
% Lorentz-Lorenz for the index of a mixture
% n1,2,3,4          Index of constituent
% phi1,1,2,3,4      Volume fractions
%
% At least two constituents required.  3 and 4 are optional.
% ***********************************************************************
function n = LorentzLorenz(n1,phi1,n2,phi2,n3,phi3,n4,phi4)

k = (n1.^2-1)./(n1.^2+2).*phi1 + (n2.^2-1)./(n2.^2+2).*phi2;

if nargin > 4
    k = k + (n3.^2-1)./(n3.^2+2).*phi3;
end

if nargin > 6
    k = k + (n4.^2-1)./(n4.^2+2).*phi4;
end

n = sqrt((1+2*k)./(1-k));

end