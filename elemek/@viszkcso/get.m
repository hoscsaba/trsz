function val = get(vena,propname)

switch propname
    case 'name'
        val=vena.name;
    case 'x'
        for i=1:vena.N+1, val(i)=vena.x(i); end
    case 'p'
        for i=1:vena.N+1, val(i)=vena.p(i); end
    case 'v'
        for i=1:vena.N+1, val(i)=vena.v(i); end
    case 'v1'
        val=vena.v(1);
    case 'vN'
        val=vena.v(vena.N+1);
   case 'p1'
        val=vena.p(1);
    case 'pN'
        val=vena.p(vena.N+1);
    case 'dt'
        val=vena.dt;
    case 'all'
        val=vena.all;
    case 'hely'
        val=vena.hely; 
    case 'ido'
        val=vena.ido;
    case 'friction'
        val=vena.friction;
    case 'A'
        val=vena.A;
    case 'a'
        val=vena.a; 
    case 'ro'
        val=vena.ro;        
    case 'fp'
        val=vena.fp;
    case 'p_roav'
        val=vena.p(1)-vena.ro*vena.a*vena.v(1);
    otherwise
        error('Invalid pipe property; name,x,p,v,p1,pN,v1,vN,dt,hely,ido,friction,A');
end
