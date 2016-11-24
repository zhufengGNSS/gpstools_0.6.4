function biass=dcbs_p1c1
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite p1-c1 dcb
% [func]   : satellite p1-c1 dcb
% [argin]  : none
% [argout] : biass = satellite p1-c1 dcbs (ns)
%                biass(n,:)=[mjd,b1,b2,b3,...]; (mjd=start date, bn=prn n dcb)
% [note]   : IGS standard analysis
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
% [history]: 06/07/08  0.1  separated parameter file
%            08/12/06  0.2  add P1-C1 according to IGSMAIL-5460
%-------------------------------------------------------------------------------
biass=[
51636,... % 2000/04/02- : IGS Mail #2744
-0.223d0, -1.027d0,  0.173d0,  1.528d0, -0.650d0,...
 0.574d0, -0.987d0, -0.801d0,  0.390d0, -1.551d0,...
-0.117d0, -9.999d9,  1.755d0,  0.574d0, -0.991d0,...
-0.674d0, -0.887d0,  0.173d0,  0.233d0, -9.999d9,...
-0.280d0, -1.564d0, -0.490d0,  0.440d0,  0.807d0,...
 1.444d0, -0.023d0, -9.999d9,  0.987d0,  1.805d0,...
-0.610d0, -9.999d9, -9.999d9, -9.999d9, -9.999d9,...
-9.999d9, -9.999d9, -9.999d9, -9.999d9, -9.999d9...

51713,... % 2000/06/18- : IGS Mail #2879
-0.254d0, -0.871d0,  0.247d0,  1.695d0, -0.580d0,...
 0.697d0, -0.907d0, -0.524d0,  0.414d0, -1.548d0,...
-0.037d0, -9.999d9,  1.831d0, -9.999d9, -1.054d0,...
-0.617d0, -0.761d0,  0.277d0,  0.270d0, -0.854d0,...
-0.457d0, -1.501d0, -0.494d0,  0.600d0,  0.687d0,...
 1.544d0,  0.013d0, -9.999d9,  0.937d0,  1.868d0,...
-0.624d0, -9.999d9, -9.999d9, -9.999d9, -9.999d9,...
-9.999d9, -9.999d9, -9.999d9, -9.999d9, -9.999d9...

51923,... % 2001/01/14- : IGS Mail #3160
 0.222d0, -0.546d0,  0.042d0,  1.294d0, -0.798d0,...
 0.625d0, -0.523d0, -0.193d0,  0.048d0, -1.002d0,...
-0.329d0, -9.999d9,  1.545d0, -0.409d0, -0.755d0,...
-9.999d9, -0.522d0, -9.999d9,  0.582d0, -0.958d0,...
-0.172d0, -1.374d0, -1.018d0,  0.459d0,  0.775d0,...
 1.077d0,  0.213d0, -0.144d0,  0.611d0,  1.745d0,...
-0.496d0, -9.999d9, -9.999d9, -9.999d9, -9.999d9,...
-9.999d9, -9.999d9, -9.999d9, -9.999d9, -9.999d9...

51986,... % 2001/03/18- : IGS Mail #3220
 0.039d0, -0.469d0,  0.104d0,  1.655d0, -0.634d0,...
 0.681d0, -0.192d0, -0.309d0,  0.317d0, -1.003d0,...
-0.161d0, -9.999d9,  1.376d0, -0.391d0, -0.238d0,...
-9.999d9, -0.695d0, -0.242d0,  0.541d0, -1.114d0,...
-0.176d0, -1.664d0, -1.033d0,  0.708d0,  0.373d0,...
 1.022d0,  0.123d0, -0.323d0,  0.511d0,  1.621d0,...
-0.426d0, -9.999d9, -9.999d9, -9.999d9, -9.999d9,...
-9.999d9, -9.999d9, -9.999d9, -9.999d9, -9.999d9...

52294,... % 2002/01/20- : IGS Mail #3674
 0.377d0, -0.611d0,  0.343d0,  1.109d0, -0.738d0,...
 0.329d0, -0.557d0, -0.061d0,  0.172d0, -1.226d0,...
 0.229d0, -9.999d9,  1.519d0, -0.279d0, -0.751d0,...
-9.999d9, -0.722d0, -0.666d0, -9.999d9, -0.953d0,...
-0.088d0, -0.626d0, -1.308d0,  0.167d0,  0.791d0,...
 0.888d0,  0.367d0, -0.217d0,  0.760d0,  2.015d0,...
-0.261d0, -9.999d9, -9.999d9, -9.999d9, -9.999d9,...
-9.999d9, -9.999d9, -9.999d9, -9.999d9, -9.999d9...

52700,... % 2003/03/02- : IGS Mail #4279
-0.067d0, -0.944d0,  0.106d0,  1.508d0, -0.802d0,...
 0.645d0, -0.916d0, -0.514d0,  0.380d0, -1.480d0,...
 0.692d0, -9.999d9,  1.503d0,  0.289d0, -0.830d0,...
-0.561d0, -0.595d0,  0.084d0, -9.999d9, -1.084d0,...
-9.999d9, -1.609d0, -0.740d0,  0.347d0,  0.720d0,...
 1.223d0, -0.023d0, -0.113d0,  0.867d0,  2.211d0,...
-0.296d0, -9.999d9, -9.999d9, -9.999d9, -9.999d9,...
-9.999d9, -9.999d9, -9.999d9, -9.999d9, -9.999d9...

52777,... % 2003/05/18- : IGS Mail #4366
-0.107d0, -1.062d0,  0.149d0,  1.535d0, -0.890d0,...
 0.596d0, -0.618d0, -0.513d0,  0.320d0, -1.658d0,...
 0.755d0, -9.999d9,  1.559d0,  0.427d0, -0.947d0,...
-0.285d0, -0.858d0,  0.085d0, -9.999d9, -1.019d0,...
-0.267d0, -1.422d0, -0.745d0,  0.345d0,  0.611d0,...
 1.322d0,  0.005d0, -0.096d0,  0.956d0,  2.065d0,...
-0.244d0, -9.999d9, -9.999d9, -9.999d9, -9.999d9,...
-9.999d9, -9.999d9, -9.999d9, -9.999d9, -9.999d9...

53039,... % 2004/02/04- : IGS Mail #4825
-0.052d0, -1.096d0,  0.015d0,  1.383d0, -0.821d0,...
 0.607d0, -0.942d0, -0.603d0,  0.392d0, -1.400d0,...
 0.487d0, -9.999d9,  1.435d0,  0.180d0, -0.926d0,...
-0.517d0, -0.811d0, -0.066d0, -9.999d9, -1.109d0,...
-0.437d0,  0.374d0, -0.426d0,  0.337d0,  0.569d0,...
 1.293d0, -0.062d0, -0.276d0,  0.785d0,  2.019d0,...
-0.333d0, -9.999d9, -9.999d9, -9.999d9, -9.999d9,...
-9.999d9, -9.999d9, -9.999d9, -9.999d9, -9.999d9...

53141,... % 2004/05/16- : IGS Mail #4937
-0.076d0, -9.999d9,  0.014d0,  1.473d0, -0.872d0,...
 0.565d0, -0.809d0, -0.563d0,  0.327d0, -1.579d0,...
 0.609d0, -9.999d9,  1.600d0,  0.311d0, -1.040d0,...
-0.460d0, -0.915d0,  0.014d0, -2.410d0, -1.084d0,...
-0.298d0,  0.652d0, -9.999d9,  0.258d0,  0.552d0,...
 1.247d0, -0.020d0, -0.162d0,  0.949d0,  2.082d0,...
-0.366d0, -9.999d9, -9.999d9, -9.999d9, -9.999d9,...
-9.999d9, -9.999d9, -9.999d9, -9.999d9, -9.999d9...

53232,... % 2004/08/15- : IGS Mail #4987
-0.053d0, -9.999d9, -0.011d0,  1.450d0, -0.944d0,...
 0.538d0, -1.228d0, -0.243d0,  0.401d0, -1.582d0,...
 0.553d0, -9.999d9,  1.621d0,  0.441d0, -0.948d0,...
-0.410d0, -0.835d0,  0.063d0, -2.323d0, -1.074d0,...
-0.345d0,  0.605d0, -0.225d0,  0.289d0,  0.570d0,...
 1.277d0,  0.001d0, -0.204d0,  0.969d0,  2.095d0,...
-0.449d0, -9.999d9, -9.999d9, -9.999d9, -9.999d9,...
-9.999d9, -9.999d9, -9.999d9, -9.999d9, -9.999d9...

53386,... % 2005/01/16- : IGS Mail #5078
-0.028d0, -0.061d0,  0.077d0,  1.334d0, -0.929d0,...
 0.664d0, -0.912d0, -0.335d0,  0.529d0, -1.567d0,...
 0.535d0, -9.999d9,  1.541d0,  0.335d0, -1.057d0,...
-0.419d0, -0.906d0,  0.098d0, -2.269d0, -1.105d0,...
-0.346d0,  0.579d0, -0.221d0,  0.152d0,  0.735d0,...
 1.247d0, -0.018d0, -0.205d0,  0.842d0,  2.017d0,...
-0.307d0, -9.999d9, -9.999d9, -9.999d9, -9.999d9,...
-9.999d9, -9.999d9, -9.999d9, -9.999d9, -9.999d9...

53646,... % 2005/10/03- : IGS Mail #5260
-0.233d0, -0.033d0, -0.163d0,  1.234d0, -0.944d0,...
 0.444d0, -1.132d0, -0.416d0,  0.366d0, -1.651d0,...
 0.591d0, -9.999d9,  1.529d0,  0.108d0, -1.326d0,...
-0.525d0,  1.391d0, -0.012d0, -2.026d0, -1.228d0,...
-0.392d0,  0.532d0,  0.130d0, -0.069d0,  0.575d0,...
 1.083d0, -0.189d0, -0.232d0,  0.643d0,  2.016d0,...
-0.073d0, -9.999d9, -9.999d9, -9.999d9, -9.999d9,...
-9.999d9, -9.999d9, -9.999d9, -9.999d9, -9.999d9...
];
