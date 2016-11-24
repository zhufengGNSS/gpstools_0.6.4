function t_ecsftoecef
sec2rad=pi/180/3600;

erp_value0=[0.5*sec2rad,-0.03*sec2rad,0.5,0.3,0.4];
tm=caltomjd([2004,1,2,3,4,5]);

[U,P,N,gmst,dUdxp,dUdyp,dUddt]=ecsftoecef(tm,erp_value0);

erp_value=erp_value0+[0.00001,0,0,0,0];
(ecsftoecef(tm,erp_value)-U)./0.00001
dUdxp

erp_value=erp_value0+[0,0.00001,0,0,0];
(ecsftoecef(tm,erp_value)-U)./0.00001
dUdyp

erp_value=erp_value0+[0,0,0.00001,0,0];
(ecsftoecef(tm,erp_value)-U)./0.00001
dUddt

