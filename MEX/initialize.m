currentpath = which('CEA.m');
[pathstr,~,~] = fileparts(currentpath);
test = CEA;
ioinp = test.input.rocket( 3 , 350, 'psia', 4.84 , 23.8 , 'paraffin' , 100, 298.15 , 'N2O', 100, 298.15);
% mex cea2.f
% data = cea2(ioinp,'\wrapper',pathstr)
