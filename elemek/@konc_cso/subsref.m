function out = subsref(buko,index)

switch index.type
 case '.'
  switch index.subs
   case 'Cd'
    out = buko.Cd;
   case 'B'
    out = buko.B;    
   case 'h0'
    out = buko.h0;
   case 'kitevo'
    out = buko.kitevo;
    
   otherwise    
    out=subsref(buko.tranziens_agelem_2csp,index);        
  end
 otherwise
  error('HIBA!!! mezohivatkozas .-al (pl. elem.ro)')
end