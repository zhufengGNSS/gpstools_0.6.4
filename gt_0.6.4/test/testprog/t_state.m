function t_state

Re=6378136.3;
x0=eletostate([Re+650000,0.001,51,0,0,0]);
t0=caltomjd([2000,1,1,0,0,0]);

%disp('state_kepler,state_2body(t=1)');
%for n=1:1000
%    [x,phi]=state_kepler(t0,1,x0,1);
%    [x,phi]=state_2body(t0,1,x0,1);
%end

disp('state_precorbit(nmax=1)-state_kepler(t=60)')
[xp,phip]=state_precorbit(t0,60,x0,1,[1,0,0]);
[xk,phik]=state_kepler(t0,60,x0,1);
[xp-xk,(xp-xk)./xk], (phip-phik)./phik

disp('state_precorbit(nmax=2)-state_geoj2(t=60)')
[xj,phij]=state_geoj2(t0,60,x0,1);
[xp,phip]=state_precorbit(t0,60,x0,1,[2,0,0]);
[xp-xj,(xp-xj)./xj], (phip-phij)./phij

%disp('state_2body-state_kepler(t=3600)')
%[x2,phi2]=state_2body(t0,3600,x0,1);
%[xk,phik]=state_kepler(t0,3600,x0,1);
%[x2-xk,(x2-xk)./xk], (phi2-phik)./phik

disp('state_precorbit(nmax=1)-state_kepler(t=3600)')
[xp,phip]=state_precorbit(t0,3600,x0,1,[1,0,0]);
[xk,phik]=state_kepler(t0,3600,x0,1);
[xp-xk,(xp-xk)./xk], (phip-phik)./phik

disp('state_precorbit(nmax=2)-state_geoj2(t=3600)')
[xj,phij]=state_geoj2(t0,3600,x0,1);
[xp,phip]=state_precorbit(t0,3600,x0,1,[2,0,0]);
[xp-xj,(xp-xj)./xj], (phip-phij)./phij

%disp('state_2body-state_kepler(t=86400)')
%[x2,phi2]=state_2body(t0,86400,x0,1);
%[xk,phik]=state_kepler(t0,86400,x0,1);
%[x2-xk,(x2-xk)./xk], (phi2-phik)./phik

%disp('state_precorbit(nmax=1)-state_kepler(t=86400)')
%[xp,phip]=state_precorbit(t0,86400,x0,1,[1,0,0]);
%[xk,phik]=state_kepler(t0,86400,x0,1);
%[xp-xk,(xp-xk)./xk], (phip-phik)./phik

%disp('state_precorbit(nmax=2)-state_geoj2(t=86400)')
%[xj,phij]=state_geoj2(t0,86400,x0,1)
%[xp,phip]=state_precorbit(t0,86400,x0,1,[2,0,0]);
%[xp-xj,(xp-xj)./xj], (phip-phij)./phij

% ‘¬“xŒv‘ª
disp('state_precorbit(t=1,10,60,300,600,3600,86400)');
for n=1:100, [x,phi]=state_precorbit(t0,1,x0,1,[8,1,0]); end
for n=1:10,  [x,phi]=state_precorbit(t0,10,x0,1,[8,1,0]); end
for n=1:10,  [x,phi]=state_precorbit(t0,60,x0,1,[8,1,0]); end
for n=1:10,  [x,phi]=state_precorbit(t0,300,x0,1,[8,1,0]); end
[x,phi]=state_precorbit(t0,600,x0,1,[8,1,0]);
[x,phi]=state_precorbit(t0,3600,x0,1,[8,1,0]);
[x,phi]=state_precorbit(t0,86400,x0,1,[8,1,0]);

