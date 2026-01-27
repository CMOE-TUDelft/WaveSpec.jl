clc
clear
close all

spec = jonswap([],[0.4 2.5],1)

figure(2)
dw = spec.w(2)-spec.w(1);
a = sqrt( spec.S*dw*2 );
plot(spec.w, a)
dw = spec.w(2) - spec.w(1)
4*sqrt(simpson(spec.w, spec.S))