function info(fojt,flag,varargin)

% flag=0 -> semmi
% flag=1 -> képernyõ
% flag=2 -> file, ami a varargin{1}-ben van
% flag=3 -> futás közbeni információ

info(fojt.tranziens_agelem_2csp,flag,varargin);

switch flag
    case 1
        fprintf('\n   %s:',class(fojt));
        fprintf('\n\tK0:              %g',fojt.K0);
        fprintf('\n\tK1:              %g',fojt.K1);        
        
    case 2
        fp=fopen(char(varargin{1}),'a');
        fprintf(fp,'\n   %s:',class(fojt));
        fprintf(fp,'\n\tK0:              %g',fojt.K0);
        fprintf(fp,'\n\tK1:              %g',fojt.K1);
        fprintf(fp,'\n-----------------------------------------------------------------\n');
        fclose(fp);
                
    otherwise
        fprintf('\n');
end