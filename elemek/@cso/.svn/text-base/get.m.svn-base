function val = get(cso,propname)

switch propname
    case 'name'
        val=cso.name;
    case 'x'
        for i=1:cso.N+1, val(i)=cso.x(i); end
    case 'p'
        for i=1:cso.N+1, val(i)=cso.p(i); end
    case 'v'
        for i=1:cso.N+1, val(i)=cso.v(i); end
    case 'v1'
        val=cso.v(1);
    case 'vN'
        val=cso.v(cso.N+1);
   case 'p1'
        val=cso.p(1);
    case 'pN'
        val=cso.p(cso.N+1);
    case 'dt'
        val=cso.dt;
    case 'hely'
        val=cso.hely; 
    case 'ido'
        val=cso.ido;
    case 'friction'
        val=cso.friction;
    case 'A'
        val=cso.A;
    case 'a'
        val=cso.a; 
    case 'ro'
        val=cso.ro;        
    case 'fp'
        val=cso.fp;
    case 'p_roav'
        val=cso.p(1)-cso.ro*cso.a*cso.v(1);
    otherwise
        error('Invalid pipe property; name,x,p,v,p1,pN,v1,vN,dt,hely,ido,friction,A');
end
