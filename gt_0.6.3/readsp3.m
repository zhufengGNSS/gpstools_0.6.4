%*------------------------------------------------------------------------------
% [system] : GpsTools
% [system] : GpsTools
% [module] : read sp3 ephemeris/clock file
% [func]   : read sp3 ephemeris/clock file
% [argin]  : file  = file name
% [argout] : epoch = start epoch [year,month,day,hour,min,sec]
%            time  = time vector (sec)
%            eph   = satellite ephemerides
%                    [x,y,z,cb(,vx,vy,vz,cd);...]
%                    x,y,z    : satellite position (ecef) (m)
%                    cb       : satellite clock bias (sec)
%                    vx,py,pz : satellite velocity (ecef) (m/sec)
%                    cd       : satellite clock drift (sec/sec)
%            sats  = satellite list
%                    ('GPSnn':GPS,'GLOnn':GLONASS,'LEOnn':LEO,'GALnn':GALILEO)
%           (accs) = accuracy of satellite orbits
%           (std)  = standard deviations for sp3c
%           (type) = file type
%                    ('':UNKNOWN,'G':GPS,'M':MIXED,'R':GLONASS,L:'LEO',E:'GALILEO)
%           (tsys) = time system
%                    ('':UNKNOWN,'GPS':GPS TIME,'UTC':UTC)
% [note]   :
% [version]: $Revision: 2 $ $Date: 06/07/08 14:21 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/03/21  0.1  mfile->mex
%            05/05/21  0.2  support argouts accs,std,type,tsys and sp3c
%-----------------------------------------------------------------------------*/

% (mex function)

