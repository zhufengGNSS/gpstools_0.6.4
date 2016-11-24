function yinterp=ntrp113(tinterp,t,y,tnew,ynew,klast,phi,psi)

yinterp=feval('C:\MATLABR11\toolbox\matlab\funfun\ntrp113',...
              tinterp,t,y,tnew,ynew,klast,phi,psi);
yinterpm=ntrp113m(tinterp,t,y,tnew,ynew,klast,phi,psi);

dv=max(max((yinterp-yinterpm)./yinterp));
if abs(dv)>eps*10, t,dv, warning('ntrp113 warning'), end
