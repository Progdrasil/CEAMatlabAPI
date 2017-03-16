% currentpath = which('CEA.m');
% [pathstr,~,~] = fileparts(currentpath);
test = CEA;
% inp = test.input.rocket( 3 , 350, 'psia', 4.84 , 23.8 , 'paraffin' , 100, 298.15 , 'N2O', 100, 298.15);
% inp = test.input.rocket( 3 , 350 ,'psia', 4.84 , 23.8 , {'paraffin' 'CH4' 'RP-1'} , [50 25 25], [298.15 298.15 298.15], {'N2O' 'O2(L)'}, [75 25],[298.15 90.1]);
% mex cea2.f
% data = cea2(test.ioinp,'\wrapper',pathstr)
test.OF = 3;
test.pressure = 350;
test.presUnit = 'psia';
test.supar = 4.84;
test.PcPe = 23.8;
test.setFuel('paraffin' , 100, 298.15);
% test.setFuel({'paraffin' 'CH4' 'RP-1'} , [50 25 25], [298.15 298.15 298.15]);
test.setOxid('N2O', 100, 298.15);
% test.setOxid({'N2O' 'O2(L)'}, [75 25],[298.15 90.1]);

ioinp = test.input.rocket;
data = test.run;