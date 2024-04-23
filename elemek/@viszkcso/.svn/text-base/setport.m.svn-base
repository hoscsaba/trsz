function out=setport(flag,vena,varargin)

N = venav.N; ro = venav.ro; a = venav.a; dzdx=0;

switch flag
    case -1
        out{1}=1;
        out{2}=1;
        
    case 0 % stacioner számítás
        error('Rugalmas csövet nem lehet stacionerben számolni!');
        
    case 1 % tranzeins számítás
        error('Nem jók a beállítások! (setport.m)');
        i   = varargin{1};
        x   = varargin{2};
        csp = varargin{3};
        dt  = varargin{4};

        kp =-(ro*9.81*a*dzdx + 0.02/2/venav.D*ro*a*venav.v(N)*abs(venav.v(N)));
        km =  ro*9.81*a*dzdx + 0.02/2/venav.D*ro*a*venav.v(2)*abs(venav.v(2));
        br = venav.p(1)   - ro*a*venav.v(1);
        ar = venav.p(N+1) + ro*a*venav.v(N+1);
        
        out{1}={-a/venav.A, 0, 0, -dt*km-br, [venav.csp1,1e5]}; 
        out{2}={ a/venav.A, 0, 0, -dt*kp-ar, [venav.csp2,1e5]};
        
    otherwise
        error('rugalmas csovet csak instac opcioval lehet hivni!');
end