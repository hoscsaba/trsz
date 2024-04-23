function val = get(csatorna,propname)

switch propname
    case 'name'
        val=csatorna.name;
    case 'x'
        for i=1:csatorna.N+1, val(i)=csatorna.x(i); end
    case 'p'
        for i=1:csatorna.N+1, val(i)=csatorna.p(i); end
    case 'v'
        for i=1:csatorna.N+1, val(i)=csatorna.v(i); end
    case 'v1'
        val=csatorna.v(1);
    case 'vN'
        val=csatorna.v(csatorna.N+1);
   case 'p1'
        val=csatorna.p(1);
    case 'pN'
        val=csatorna.p(csatorna.N+1);
    case 'dt'
        val=csatorna.dt;
    case 'hely'
        val=csatorna.hely; 
    case 'ido'
        val=csatorna.ido;
    case 'friction'
        val=csatorna.friction;
    case 'A'
        val=csatorna.A;
    case 'a'
        val=csatorna.a; 
    case 'ro'
        val=csatorna.ro;        
    case 'fp'
        val=csatorna.fp;
    case 'p_roav'
        val=csatorna.p(1)-csatorna.ro*csatorna.a*csatorna.v(1);
    otherwise
        error('Invalid pipe property; name,x,p,v,p1,pN,v1,vN,dt,hely,ido,friction,A');
end
