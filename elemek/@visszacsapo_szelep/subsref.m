function out = subsref(sziv,index)

switch index.type
 case '.'
  switch index.subs
   otherwise
    out=subsref(sziv.tranziens_agelem_2csp,index);      
  end      
 otherwise 
  error('HIBA!!! mezohivatkozas .-al (pl. elem.ro)')
end