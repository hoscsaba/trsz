function mar = set(mar,varargin)

property_argin = varargin;
while length(property_argin) >= 3,
    i    = property_argin{1};
    prop = property_argin{2};
    val  = property_argin{3};
    property_argin = property_argin(4:end);
    
    switch prop
        case 'Q'
            temp = mar.elemek{i};
            temp.Q = val;
            mar.elemek{i} = temp;
            
        case 'p'
            ize=1;
            while ~(mar.csp{ize}{1}==i)
                if ize==length(mar.csp), error(['Nincs ilyen szamu csomopont: ',num2str(hol),' a ',mar.nev,' nevu merev alrendszerben']), end
                ize=ize+1;
            end
            mar.csp{ize}{4}=val;
            
        case 'resfile'
            temp=mar.elemek{i};
            temp.resfile=val;
            mar.elemek{i} = temp;
            
        case 'save_level'
            mar.save_level=val;
            
        case 'plot_iter'
            mar.plot_iter=val;
        otherwise
            error('Invalid property')
    end
end