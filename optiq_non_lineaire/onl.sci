// ****************************************************************
//  Bibliothèque sommaire associée au cours d'optique non-linéaire
//  (c) 2014 Ecole Polytechnique
//  Version 1.3
// ****************************************************************

// Modification des limites de l'axe des abscisses d'un graphe
// Doit être appelé juste après plot
function xlim(x)
  ax = gca();
  ax.data_bounds(1,1) = x(1);
  ax.data_bounds(2,1) = x(2);
endfunction;

// Modification des limites de l'axe des ordonnées d'un graphe
// Doit être appelé juste après plot
function ylim(y)
  ax = gca();
  ax.data_bounds(1,2) = y(1);
  ax.data_bounds(2,2) = y(2);
endfunction;

// Phase (en radians) d'un nombre complexe
function result = angle(z)
  z = z + bool2s(z==0); // Pour éviter une erreur si z=0.
  result = atan(imag(z),real(z))
endfunction;

// Supprime par continuité les sauts de 2pi dans la phase
function result = unwrap(phi)
  result = phi - [0 2*%pi*cumsum(floor(.5+diff(phi)/2/%pi))];
endfunction;

// Prend l'inverse terme à terme d'un vecteur
// Dans Scilab, si x est une matrice n x 1, 1./x donne une matrice
// 1 x n telle que le produit de la matrice 1xn par la matrice
// initiale nx1 donne 1.
// Il faudrait en principe utiliser 1 ./ x avec un espace entre le
// 1 et le .
// inverse(x) paraît plus lisible que 1 ./x car on a vite fait
// d'oublier l'espace.
function result = inverse(x)
  result = 1 ./x;
endfunction;

// Prend le carré d'une grandeur (plus rapide que .^2)
function result = sqr(x)
  result = x .* x;
endfunction;

// Fonction de Heaviside
function result = heaviside(x)
  // (x>0) est un tableau de boolean
  // transformé en 0 ou 1 à l'aide de bool2s
  result = bool2s(x>0);
endfunction;

// Ecart quadratique d'une grandeur x pondérée par l'amplitude de probabilité field
function result = ecartQuad(x,field)
  n2 = sum(sqr(abs(field)));
  xMoy = sum(x.*sqr(abs(field)))/n2;
  x2Moy = sum(x.*x.*sqr(abs(field)))/n2;
  result = sqrt(x2Moy-sqr(xMoy));
endfunction;

// Choix de l'échelle de couleur utilisée pour les représentations en fausses couleurs
// 0 pour niveaux de gris, 1 pour jetcolormap
function setColorMap(i)
  hcf = gcf();
  select i
    case 0 then
      hcf.color_map = graycolormap(256);
    case 1 then
      hcf.color_map = jetcolormap(256);
  end;
endfunction;

// Lecture d'un fichier ascii au format csv (comma separated values)
function result = csvread(fileName)
  fd = mopen(fileName);
  str = mgetl(fd);
  [n dummy] = size(str); // n est le nombre de lignes
  [m dummy] = size(tokens(str(1),",")); // m est le nombre de colonnes
  for i=1:n
    line = tokens(str(i),",");
    for j=1:m
      result(i,j) = sscanf(line(j),"%e");
    end;
  end;
  mclose(fd);
endfunction;

// Décale un tableau array de la quantité spécifiée shift
// en répétant le premier (resp. dernier) point si le décalage est
// positif (resp. négatif)
function result = shiftArray(array, shift)
    result = ones(array);
    if (shift>0) then
        result = result.*array(1);
        result(shift+1:$)=array(1:$-shift);
    elseif (shift<0) then
        result = result.*array($);
        result(1:$+shift)=array(1-shift:$);
    else // shift = 0
        result = array;
    end
endfunction

// Détermine les maxima des données spécifiées
// L'argument optionnel width spécifie la largeur (en pixels) de la zone
// dont les points trouvés doivent être maximum
// Par défaut width = 10, ce qui signifie que chaque maximum est le point
// le plus haut d'une zone de 21 pixels de large
function result = findMaxima(data,width)
    if (argn(2)==1) then
      width = 10; // Valeur par défaut du deuxième argument optionnel
    end
    threshold = zeros(data);
    for i=1:width
      threshold = max(threshold,shiftArray(data,i),shiftArray(data,-i));
    end
    result = find(data>threshold);
endfunction

// Transformée de Fourier avec rotation pour centrer la fréquence nulle
function result = ft(data)
  result = fftshift(fft(fftshift(data)));
endfunction;

// Transformée de Fourier inverse avec rotation pour centrer la fréquence nulle
function result = ift(data)
  result = fftshift(ifft(fftshift(data)));
endfunction;

// ftAxis
// Création de deux tableaux (fréquence et temps) calibrés l'un par rapport à l'autre pour
// une FFT. L'espacement entre deux points de fréquence est l'inverse de la largeur totale
// dans l'espace des temps. Inversement, l'espacement entre deux points de temps est l'inverse
// de la largeur totale en fréquence.
// Les axes de fréquence et de temps sont supposés centré sur zéro
//
// nPoints : Nombre de points (un nombre obligatoirement pair, de préférence une puissance de 2)
// nuMax   : Valeur maximale de la fréquence (ie moitié de la largeur totale)
// nu : Axe des fréquences. Le point d'indice nPoints/2 vaut toujours zéro
// t  : Axe des temps. Le point d'indice nPoints/2 vaut toujours zéro
function [nu, t] = ftAxis(nPoints, nuMax)
  deltaNu = 2*nuMax/nPoints;
  deltaT = 1/(2*nuMax);
  nu = -nuMax:2*nuMax/nPoints:nuMax-(2*nuMax/nPoints);
  t = -nPoints/2*deltaT:deltaT:(nPoints/2-1)*deltaT;
endfunction;

// Formule de Sellmeier pour le calcul de l'indice de réfraction
// n^2 = 1 + B1/(1-C1/lambda^2) + B2/(1-C2/lambda^2) + B3/(1-C3/lambda^2)
// lambda s'exprime en microns
function result = sellmeier(lambda,b1,b2,b3,c1,c2,c3)
  lm2 = lambda.^-2;
  result = sqrt(1+b1./(1-c1*lm2)+b2./(1-c2*lm2)+b3./(1-c3*lm2));
endfunction;

// Dérivée de la formule de Sellmeier pour le calcul de l'indice de groupe
// lambda s'expriem en microns
function result = sellmeierGroupe(lambda,b1,b2,b3,c1,c2,c3)
    lm2 = lambda.^(-2);
    lm3 = lambda.^(-3);
    ndn = -b1*c1.*lm3./sqr(1-c1*lm2)-b2*c2.*lm3./sqr(1-c2*lm2)-b3*c3.*lm3./sqr(1-c3*lm2);
    n = sellmeier(lambda,b1,b2,b3,c1,c2,c3);
    result = n-ndn./n.*lambda;
endfunction

// Indices ordinaire et extraordinaire de quelques matériaux courants
// La longueur d'onde lambda s'exprime en microns
function [no, ne] = indice(cristal,lambda)
  lambda2 = sqr(lambda);
  select cristal
    case 'BBO' then
      no = sqrt(2.7405+.0184./(lambda2-.0179)-.0155 *lambda2);
      ne = sqrt(2.3730+.0128./(lambda2-.0156)-.0044 *lambda2);
      case 'AgGaS2' then
      no = sqrt(3.40684+2.40065*lambda2./(lambda2-.09311)+2.06248*lambda2./(lambda2-950));
      ne = sqrt(3.60728+1.94792*lambda2./(lambda2-.11066)+2.24544*lambda2./(lambda2-1030.7));
    case 'KDP' then
      no = sqrt(2.259276 + 0.01008956./(lambda2 - 0.012942625) ...
                         + 13.005522*lambda2./(lambda2 - 400));
      ne = sqrt(2.132668 + 0.008637494./(lambda2 - 0.012281043) ...
                         + 3.2279924*lambda2./(lambda2 - 400));
    case 'LiNbO3' then
      no = sqrt(1 + 2.6734*lambda2./(lambda2 - 0.01764)+ 1.2290*lambda2./(lambda2 - 0.05914) + 12.614*lambda2./(lambda2 - 474.60));
      ne = sqrt(1 + 2.9804*lambda2./(lambda2 - 0.02047)+ 0.5981*lambda2./(lambda2 - 0.0666) + 8.9543*lambda2./(lambda2 - 416.08));
    case 'LiIO3' then
      no = sqrt(2.03132 + 1.37623*lambda2./(lambda2 - 0.0350823)+ 1.06745*lambda2./(lambda2 - 169));
      ne = sqrt(1.83086 + 1.08807*lambda2./(lambda2 - 0.0313810)+ 0.0554582*lambda2./(lambda2 - 158.76));
    case 'SiO2' then
      // http://cvilaser.com/Common/PDFs/Dispersion_Equations.pdf
      no = sellmeier(lambda,.6961663,.4079426,.8974794,4.67914826E-3,1.35120631E-2,9.79340025E1);
      ne = no;
    case 'SF10' then
      no = sellmeier(lambda,1.61625977,0.259229334,1.07762317,0.0127534559,0.0581983954,116.60768);
      ne = no;
    case 'BK7' then
      no = sellmeier(lambda,1.03961212,0.231792344,1.01046945,0.00600069867,0.0200179144,103.560653);
      ne = no;
    case 'CaF2' then
      // http://www.us.schott.com/lithotec/english/download/caf2_june_2006_final_us.pdf
      no = sellmeier(lambda,.6188140,.4198937,3.426299,2.759866E-3,1.061251E-2,1.068123E3);
//      no = sellmeier(lambda,.567588800,.471091400,3.84847230,.00252642999,.0100783328,1200.55597);
      ne = no;
    else
      no = 1; ne = 1;
  end;
endfunction;

// Indice extraordinaire effectif pour un angle theta donné
function result = neTheta(cristal, theta, lambda)
  [no, ne] = indice(cristal, lambda);
  result = (sqr(cos(theta)./no)+sqr(sin(theta)./ne)).^(-.5);
endfunction;

// Indice de groupe de quelques matériaux courants
function result = indiceGroupe(cristal,lambda)
  lambda2 = sqr(lambda);
  select cristal
    case 'SiO2' then
      // http://cvilaser.com/Common/PDFs/Dispersion_Equations.pdf
      result = sellmeierGroupe(lambda,.6961663,.4079426,.8974794,4.67914826E-3,1.35120631E-2,9.79340025E1);
    case 'SF10' then
      result = sellmeierGroupe(lambda,1.61625977,0.259229334,1.07762317,0.0127534559,0.0581983954,116.60768);
    case 'BK7' then
      result = sellmeierGroupe(lambda,1.03961212,0.231792344,1.01046945,0.00600069867,0.0200179144,103.560653);
    case 'CaF2' then
      // http://www.us.schott.com/lithotec/english/download/caf2_june_2006_final_us.pdf
      result = sellmeierGroupe(lambda,.6188140,.4198937,3.426299,2.759866E-3,1.061251E-2,1.068123E3);
//      no = sellmeier(lambda,.567588800,.471091400,3.84847230,.00252642999,.0100783328,1200.55597);
    else
      result = 1;
  end;
endfunction;

// Initialisation des paramètres par défaut des représentations graphiques
hda = gda();               // Handle du repère par défaut (default axes)
hda.font_size = 5;         // Augmente la taille de police des axes
hda.title.font_size = 5;   // Augmente la taille du titre
hda.x_label.font_size = 5; // Axe des x
hda.y_label.font_size = 5; // Axe des y
// Pour projection ou figure article, utiliser paramètres ci-dessous
hda.font_size = 5;         // Augmente la taille de police des axes
hda.thickness = 2;         // Epaisseur par défaut des lignes
