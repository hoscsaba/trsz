function dt = update_dt(csatorna)

if strcmp(csatorna.dttype,'auto')
    dt=zeros(size(csatorna.y));
    for i=1:length(csatorna.y)
        AA = get_A(csatorna,csatorna.y(i));
        BB = get_B(csatorna,csatorna.y(i));
        hull_seb = sqrt(9.81*AA/BB);
        dtv(i) = csatorna.L/csatorna.N/(hull_seb+abs(csatorna.v(i)));
    end
    dt=min(2.0*min(dtv),5);
elseif ~isempty(str2double(csatorna.dttype))
    dt=str2double(csatorna.dttype);
else
    error(['A ',csatorna.tranziens_agelem_2csp.nev,' csatorna dttype valtozoja: ',csatorna.dttype,'. Megengedett ertek: auto vagy egy szam']);
end
