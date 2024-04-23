function out=setport(flag,cso,varargin)

N = csov.N; ro = csov.ro; a = csov.a; dzdx=0;

switch flag
    case -1
        out{1}=1;
        out{2}=1;
        
    case 0 % stacioner számítás
        error('Rugalmas csövet nem lehet stacionerben számolni!');
        
    case 1 % tranzeins számítás
        i   = varargin{1};
        x   = varargin{2};
        csp = varargin{3};
        dt  = varargin{4};

        kp =-(ro*9.81*a*dzdx + 0.02/2/csov.D*ro*a*csov.v(N)*abs(csov.v(N)));
        km =  ro*9.81*a*dzdx + 0.02/2/csov.D*ro*a*csov.v(2)*abs(csov.v(2));
        br = csov.p(1)   - ro*a*csov.v(1);
        ar = csov.p(N+1) + ro*a*csov.v(N+1);
        
        out{1}={-a/csov.A, 0, 0, -dt*km-br, [csov.csp1,1e5]}; 
        out{2}={ a/csov.A, 0, 0, -dt*kp-ar, [csov.csp2,1e5]};
        
    otherwise
        error('rugalmas csovet csak instac opcioval lehet hivni!');
end