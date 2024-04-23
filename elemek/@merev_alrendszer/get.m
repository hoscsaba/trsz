function val = get(mar,mit,hol)

switch mit
    case 'p'
        temp = mar.csp;
        ize = 1;
        while ~(temp{ize}{1}==hol)
            if ize == length(mar.csp)
                error(['Nincs ilyen csomopont: ',hol,' a ',mar.nev,' nevu merev alrendszerben'])
            end
            ize = ize+1;
        end
        val = temp{ize}{4}+temp{ize}{2}*9.81*1000*0;
        
    case 'pnev'
        temp = mar.csp;
        ize = 1;
        while ~strcmp(temp{ize}{6},hol)
            if ize == length(mar.csp)
                error(['Nincs ilyen csomopont: ',hol,' a ',mar.nev,' nevu merev alrendszerben'])
            end
            ize = ize+1;
        end
        val = temp{ize}{4}+temp{ize}{2}*9.81*1000*0;
        
        %% Megadott nyomasu pont, szerintem nem kell.
        %temp2 = mar.elemek;
        %if isa(temp2{ize},'nyomas')
        %    val = temp2{ize}.p
        %end
        %fprintf('A helyes nyomas pedig: %8.3f\n',temp2{ize}.p);
        
    case 'Q'
        temp = mar.elemek;
        val = temp{hol}.Q;
        
    case 'm'
        temp=mar.elemek;
        val=temp{hol}.Q*temp{hol}.ro;
        
    case 'dp'
        pp=mar.csp;
        temp=mar.elemek;
        temp2=temp{hol};
        csp=temp2.csp;
        if length(csp)==2
            cspe=csp(1);cspv=csp(2);
            val=pp{cspe}{4}-pp{cspv}{4};
        else
            val=0;
        end
        
        %  val=temp{hol}.csp*temp{hol}.ro;
        
    otherwise
        error('Invalid property')
end
