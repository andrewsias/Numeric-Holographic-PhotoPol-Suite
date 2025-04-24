%==========================================================================
% CallHoloPolymer.  Calculate and render results of HoloPolymer function
%==========================================================================

nano  = 10^-9;  micro = 10^-6;  milli = 10^-3;

Rotation = 0.0;

[in,out] = NumericalHoloPolymer('OutSaveFile',"",'InSaveFile',"",'CalcBragg',true,...
    'Z',50*micro,'ThetaIncWrt',0.1346+Rotation,'ThetaDifWrt',-0.1346+Rotation,...
    'lambda0Red',660*nano,'w0red',40*micro,'Nz',100,'NBraggAngles',101,'Distort',false ,'deltan',.001);

if out.Writing
    out = TheoryHoloPolymer(in,out);
end

RenderHoloPolymer(in,out);

