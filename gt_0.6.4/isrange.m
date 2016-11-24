%--------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : range between satellites
% [func]   : calculate range between two satellites
% [argin]  : state1 = satellite no.1 position/velocity(m) (eci)
%            state2 = satellite no.2 position/velocity(m) (eci)
%            dts1   = satellite no.1 clock bias(sec)
%            dts2   = satellite no.2 clock bias(sec)
%            corrlightt = light-time/receiver clock bias correction flag (1:on)
% [argout] : rs1    = satellite no.1 tx position (m) (eci)
%            rs2    = satellite no.2 tx position (m) (eci)
%            range  = range between satellite no.1 and no.2
%            (drds1)= partial derivatives of range by satellite no.1 state
%            (drds2)= partial derivatives of range by satellite no.2 state
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/02/11   0.1  new
%-------------------------------------------------------------------------------

% (mex function)

