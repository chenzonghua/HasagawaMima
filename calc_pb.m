%    ≤¥À…¿®∫≈hcs={u,v}           %
function pb=calc_pb(u,v)
global dx dy
global im ip jp jm ic jc
pb = (3*u(ip,jc)-2*u(ic,jc)-u(im,jc))/4/dx.*(3*v(ic,jp)-2*v(ic,jc)-v(ic,jm))/4/dy ...
    -(3*v(ip,jc)-2*v(ic,jc)-v(im,jc))/4/dx.*(3*u(ic,jp)-2*u(ic,jc)-u(ic,jm))/4/dy;