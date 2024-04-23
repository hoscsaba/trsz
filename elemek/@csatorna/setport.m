function out=setport(flag,csatorna,varargin)

N = csatorna.N; 
ro = csatorna.ro; 
a = csatorna.a; 
dzdx=0;

switch flag
    case -1
        out{1}=1;
        out{2}=1;
        
    case 0 % stacioner számítás
        error('Rugalmas csatornat nem lehet stacionerben számolni!');
        
    case 1 % tranzeins számítás
        
        i   = varargin{1};
        x   = varargin{2};
        csp = varargin{3};
        dt  = varargin{4};

        kp =-(ro*9.81*a*dzdx + 0.02/2/csatorna.D*ro*a*csatorna.v(N)*abs(csatorna.v(N)));
        km =  ro*9.81*a*dzdx + 0.02/2/csatorna.D*ro*a*csatorna.v(2)*abs(csatorna.v(2));
        br = csatorna.p(1)   - ro*a*csatorna.v(1);
        ar = csatorna.p(N+1) + ro*a*csatorna.v(N+1);
        
        out{1}={-a/csatorna.A, 0, 0, -dt*km-br, [csatorna.csp1,1e5]}; 
        out{2}={ a/csatorna.A, 0, 0, -dt*kp-ar, [csatorna.csp2,1e5]};
        
    otherwise
        error('rugalmas csatornat csak instac opcioval lehet hivni!');
end