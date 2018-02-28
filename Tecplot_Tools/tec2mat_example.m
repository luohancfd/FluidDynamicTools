clear
clc
[zone1,VARlist1] = tec2mat('./tec2mat_example/bc.dat','debug');
[zone2,VARlist2] = tec2mat('./tec2mat_example/fvs.dat','debug');
[zone3,VARlist3] = tec2mat('./tec2mat_example/N2N.dat','passive','debug');
[zone4,VARlist4] = tec2mat('./tec2mat_example/streamline.dat','debug');