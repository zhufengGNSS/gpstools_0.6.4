function [x,P]=gt_udstate(td,t0,t1,x,ix,P,state,prn,prm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : temporal update of states
% [func]   : temporal update of states
% [argin]  : td      = date(mjd-gpst)
%            t0,t1   = previous time/updated time
%            x,ix,P  = states/states index/covariances
%            state   = a priori states
%            prn     = process noise standard deviations
%            prm     = processing parameters struct
% [argout] : x,P     = updated states
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/03/26  0.1  separated from gpsestd.m
%-------------------------------------------------------------------------------

% update satellite states
if ~isempty([ix.sata{:}])
    for n=1:length(prm.sats)
        i=[ix.sato{n},ix.sats{n}]; j=ix.satc{n};
        
        % satellite orbit
        if ~isempty(i)&~isnan(x(i(1)))
            [x(i),T]=orbitmodel(td,t0,t1,x(i),prm.sats{n},prm.sat.srp(n,:),...
                                state.erp,prm);
            P(i,:)=T*P(i,:); P(:,i)=P(:,i)*T';
            if state.ecls(n)==2, prn(i)=prn(i)*prm.eclnprns(n); end
        end
        % satellite clock
        if ~isempty(j)&~isnan(x(j(1)))
            switch prm.sclkmodel
            case 0, % white-noise
                x(j)=nan;
            case 1, % gauss-marcov
                T=[1,(t1-t0)/3600;0,1]; x(j)=T*x(j);
                P(j,:)=T*P(j,:); P(:,j)=P(:,j)*T';
            end
        end
    end
end
% update station states
for n=1:length(prm.rcvs)
    i=ix.rcvp{n}; j=ix.rcvc{n}; k=ix.rcvz{n};
    
    % receiver position
    if ~isempty(i)&~isnan(x(i(1)))
        switch prm.rposmodel
        case 1, % kinematic
            x(i)=nan;
        case 2, % dynamic
            T=[eye(3),eye(3)*(t1-t0);zeros(3),eye(3)];
            x(i)=T*x(i);
            P(j,:)=T*P(j,:); P(:,j)=P(:,j)*T';
        case 3, % satellite
            [x(i),T]=orbitmodel(td,t0,t1,x(i),prm.rcvs{n},[],state.erp,prm);
            P(i,:)=T*P(i,:); P(:,i)=P(:,i)*T';
        end
    end
    % receiver clock
    if ~isempty(j)&~isnan(x(j(1)))
        switch prm.rclkmodel
        case 0, % white-noise
            x(j)=nan;
        case 1, % gauss-marcov
            T=[1,(t1-t0)/3600;0,1]; x(j)=T*x(j);
            P(j,:)=T*P(j,:); P(:,j)=P(:,j)*T';
        end
    end
    % tropospheric zpd
    if ~isempty(k)&~isnan(x(k(1)))
        switch prm.zpdmodel
        case 1, % gauss-marcov
            T=[1,(t1-t0)/3600;0,1]; x(k)=T*x(k);
            P(k,:)=T*P(k,:); P(:,k)=P(:,k)*T';
        end
    end
end
% add process noises to covariances
i=find(~isnan(x')&prn>0);
if ~isempty(i), j=sub2ind(size(P),i,i); P(j)=P(j)+abs(t1-t0)*prn(i).^2; end

% satellite orbit dynamics -----------------------------------------------------
function [x,T]=orbitmodel(td,t0,t1,x,sat,srp,erp,prm)
global utc_tai erp_value
utc_tai=prm_utc_tai(td+t0/86400,1);
tu=td+(t0+19+utc_tai)/86400;
erp_value=erp;
prm.obt.satprm=srp;
prm.obt.nutmodel=prm.nutmodel;
[x,T]=state_precorbit(tu,t1-t0,x,sat,prm.obt);
