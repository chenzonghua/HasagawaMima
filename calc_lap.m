%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Computes the explicit Laplacian                %
%      Delta(u^n)                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hc=calc_lap(u)

global dx dy
global im ip jp jm ic jc

hc = (u(ip,jc)-2*u(ic,jc)+u(im,jc))/(dx*dx)...
   + (u(ic,jp)-2*u(ic,jc)+u(ic,jm))/(dy*dy);
