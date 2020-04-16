     close all; clear; clc;
     format long e;
     global dx dy Lx Ly;
     global nxm nym ;
     global ip im jp jm ic jc;
     nx=257;
     ny=257;
     dt=0.1;
     vd=1;
     Lx=8*pi; Ly=8*pi;
      nxm=nx-1  ;           nym=ny-1;
      dx=Lx/nxm ;           dy=Ly/nym;
      ic=1:nxm;             jc=1:nym; 
      xc=(ic-1)*dx-Lx/2 ;        yc=(jc-1)*dy-Ly/2 ;
      xm=(ic-0.5)*dx-Lx/2;       ym=(jc-0.5)*dy-Ly/2;
      ip=ic+1; ip(nxm)=1;   jp=jc+1; jp(nym)=1;
      im=ic-1; im(1)=nxm;   jm=jc-1; jm(1)=nym;
                             
[xx,yy]=meshgrid(xm,ym);% centers of the cells for visualization

        
%=========================================%
%           Initial condition             %
%=========================================%
phi  =0.01*exp(-0.01*(xx.^2+yy.^2));   
epl  =0.001;
phiT =cos(xx+yy);
%===============================%
% Optimization of the ADI method%
%===============================%      
      rr=1;
      bx=rr*dt/(dx*dx);
      by=rr*dt/(dy*dy);
[amix,apix,alphx,xs2x]=ADI_init(-bx*ones(1,nxm),(1+2*bx)*ones(1,nxm),-bx*ones(1,nxm));
[amiy,apiy,alphy,xs2y]=ADI_init(-by*ones(1,nym),(1+2*by)*ones(1,nym),-by*ones(1,nym));
cd('E:\matlab work 2014\Hasagawamima');%新文件夹路径
dirname='画图1';%新的文件夹名
Order=['mkdir ' dirname];%创建命令
system(Order);%创建文件夹
%===============================%
%   Time loop                   %
%===============================%
niter=0; temps=0;
nitermax=1000;   
tc=cputime; % to estimate the computational CPU time: initialization
   while (niter<=nitermax)
          niter=niter+1;
          temps=temps+dt;
          lapphi=calc_lap(phi);
          phi1=phi-dt*(vd*calc_dy(phi)+calc_pb(lapphi,phi)-epl.*phiT);
          du1 = ADI_step(amix,apix,alphx,xs2x,phi1');
	      phi = ADI_step(amiy,apiy,alphy,xs2y,du1');

figure(2)
clf;
subplot(1,2,1)
pcolor(xc,yc,phi);
shading('flat');
title(['\phi',' t=',num2str(temps)]);
xlabel('x');
ylabel('y');


subplot(1,2,2)
surfc(xc,yc,phi);
shading('flat');
title('\phi');
xlabel('x');
ylabel('y');
drawnow;
if mod(niter,20) == 0
    saveas(gcf,['E:\matlab work 2014\Hasagawamima\画图1\','psi-phi-epl-ome-',num2str(1000+niter),'.jpg']);
end
    end;
    

