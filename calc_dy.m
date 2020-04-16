function dudy=calc_dy(u)
%对y求一阶导数
global dx dy
global im ip jp jm ic 
dudy = (u(ic,jp)-u(ic,jm))/dy/2;%2阶精度
