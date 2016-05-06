// modele de polarization non-lineaire
// utilise la fonction exp(x)-1
function result = polarization(field)
    result = exp(field)-1;
endfunction

t = -1.5:.001:1.5;
field = cos(2*%pi*t);
