%==========================================================================
% Explore OCT idea. RRM 6/22
%==========================================================================

nano  = 10^-9;  micro = 10^-6;  milli = 10^-3;

%Rotation = [-.1,0,.1];
Rotation = 0.0;

[in,out] = NumericalHoloPolymer('OutSaveFile',"",'InSaveFile',"",...
    'Z',50*micro,'tau',0.25,'CalcBragg',true,...
    'lambda0Red',660*nano,'w0red',40*micro,'Nz',100,'Distort',false);

out = TheoryHoloPolymer(in,out);

RenderHoloPolymer(in,out);

% out.Ediff_fxfyZtheta