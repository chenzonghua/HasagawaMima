function dudy=calc_dy(u)
%��y��һ�׵���
global dx dy
global im ip jp jm ic 
dudy = (u(ic,jp)-u(ic,jm))/dy/2;%2�׾���
