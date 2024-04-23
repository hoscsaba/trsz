function info(cso,flag,varargin)

% flag=0 -> semmi
% flag=1 -> képernyõ
% flag=2 -> file, ami a varargin{1}-ben van
% flag=3 -> futás közbeni információ

info(cso.tranziens_agelem_2csp,flag,varargin)

switch flag
    case 1
        fprintf('\n   %s:',class(cso));
        fprintf('\n\tátmérõ:          %g [m]',cso.D);
        fprintf('\n\thossz:           %g [m]',cso.L);        
        fprintf('\n\tlambda:          %g\n',cso.lambda);        
        
    case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:',class(cso));
        fprintf(fp,'\n\tátmérõ:          %g [m]',cso.D);
        fprintf(fp,'\n\thossz:           %g [m]',cso.L);        
        fprintf(fp,'\n\tlambda:          %g',cso.lambda);    
        fprintf(fp,'\n-----------------------------------------------------------------\n');
        fclose(fp);
                
    otherwise
        fprintf('\n');
end