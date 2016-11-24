function t_state_lt

x0=eletostate([42164000,0.1,45,327,270,0]);
t0=caltomjd([2008,1,1,12,0,0]);
isat=1;

t=1;
[x1,phi1]=state_2body(t0,t,x0,isat);
[x2,phi2]=state_lt(t0,t,x0,isat);
x2
x1
phi2
phi1

function [x,phi]=state_lt(t0,t,x0,isat)
GMe=3.986004415E+14;
r=x0(1:3); v=x0(4:6); rr=norm(r);
a=-GMe/rr^3*r;
dadr=-GMe/rr^3*(eye(3)-3*r*r'/rr^2);
rr=norm(x0(1:3));
x(1:3,1)=r+v*t+a/2*t^2;
x(4:6,1)=v+a*t;
phi=[eye(3)+dadr*t^2/2,eye(3)*t;dadr*t,eye(3)];

