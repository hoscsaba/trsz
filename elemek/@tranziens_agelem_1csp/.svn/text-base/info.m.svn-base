function info(trag1,flag,varargin)

% flag=0 -> semmi
% flag=1 -> képernyõ
% flag=2 -> file, ami a varargin{1}-ben van
% flag=3 -> futás közbeni információ

switch flag
    case 1
        fprintf('\nElem: %s\n   %s',trag1.nev,class(trag1));
        fprintf('\n\tcsomópont:        %d',trag1.csp);
        fprintf('\n\tnyomás:           %g [Pa]',trag1.p);
        fprintf('\n\ttérfogatáram:     %g [m^3/s]\n',trag1.Q);
        fprintf('\n\tfolyadék:        %s',trag1.folynev);
        fprintf('\n\t  sûrûség:       %g [kg/m^3]',trag1.ro);
        fprintf('\n\t  rug.mod:       %g [Pa]',  trag1.B);
%        fprintf('\n\t  kin.viszk.:    %g [m^2/s]',trag1.nu);
%        fprintf('\n\t  din.viszk.:    %g [kgs/m^2]\n',trag1.mu);
        
    case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\nElem: %s\n   %s',trag1.nev,class(trag1));
        fprintf(fp,'\n\tcsomópont:        %d',trag1.csp);
        fprintf(fp,'\n\tnyomás:           %g [Pa]',trag1.p);
        fprintf(fp,'\n\ttérfogatáram:     %g [m^3/s]',trag1.Q);        
        fprintf(fp,'\n\tfolyadék:        %s',trag1.folynev);
        fprintf(fp,'\n\t  sûrûség:       %g [kg/m^3]',trag1.ro);
        fprintf(fp,'\n\t  rug.mod:       %g [Pa]',  trag1.B);
%        fprintf(fp,'\n\t  kin.viszk.:    %g [m^2/s]',trag1.nu);
%        fprintf(fp,'\n\t  din.viszk.:    %g [kgs/m^2]\n',trag1.mu);     
        fclose(fp);
        
    case 3
        t=varargin{1};
        fnev=char(strcat(trag1.nev,'.res'));
        if fopen(fnev,'r')==-1
            fp = fopen(fnev,'w');
            fprintf(fp,'     t [s]             p[Pa]             Q[m3/s]\n');
        else, fp = fopen(fnev,'a'); end   
        fprintf(fp,'   %5.3e       %+5.3e       %+5.3e\n',t,trag1.p,trag1.Q);
        fclose(fp);
        
    otherwise
        fprintf('\n');
end