%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : square root convergence filter
% [func]   : calculate measurement update rule by square root convergence filter
% [argin]  : x   = prefit states
%            P   = prefit covariences matrix
%            dz  = residuals
%            G   = jacobian matrix
%            R   = measurement noise corvariences matrix
%           (nsig) = outlier threshold (sigma,0:no exclusion)
% [argout] : xu  = postfit states
%            Pu  = postfit covarience matrix
%            vz  = valid residulas (1:valid,0:invalid/excluded)
% [note]   : use intel MKL CBLAS/LAPACK library
% [version]: $Revision: 1 $ $Date: 06/03/25 22:04 $
% [history]: 04/07/22   0.1  new
%-------------------------------------------------------------------------------
