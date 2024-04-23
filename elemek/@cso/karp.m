function out = karp(cso,dt,dtmin)

if dt<0
    if abs(dt)<dtmin
        dt=abs(dt);
    else
        dt
        error('Negativ dt, WTF???')
    end
end

x_interp=cso.L-cso.dx*dt/cso.dt;
p=interp1(cso.x,cso.p,x_interp); 
v=interp1(cso.x,cso.v,x_interp); 
ro = cso.tranziens_agelem_2csp.ro;

val1=p+ro*cso.a*v;
val2= - dt*( ro*cso.a*9.81*(cso.h(end)-cso.h(end-1))/cso.dx + ro*cso.a*cso.lambda/2/cso.D*v*abs(v));
val=val1+val2;

% Az elso 4 szam a fontos, a tobbi csak info
out=[1e5,  cso.a/cso.A, val,  1, val, p, v, val1, val2];

%fprintf('\n %s cs� �llapota:\n',cso.tranziens_agelem_2csp.nev);
%disp([cso.x; cso.p/1e5; cso.v]);
%fprintf('\n Ide kell C+ ment�n interpol�lni: dt=%g  cso.dt=%g  cso.dx=%g  cso.a=%g  x=%g,  p=%g,  v=%g\n',dt,cso.dt,cso.dx,cso.a,x_interp,p/1e5,v);
%pause