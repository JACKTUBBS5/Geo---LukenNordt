options nodate nonumber ps=200 ls=80;
title 'SoilGeoChem';

data geo; set sasuser.geo_MOL; run;
*proc contents; run;
data geo; set geo;
    if ph_H2O > 0;
    logSB = log_SB_;
    mAC = mAl__Al_Ca_;		
	keep ph_H2O logSB mAC;
run;
*proc contents; run;
proc sgplot data=geo;
scatter y=ph_H2O x=logSB;
loess y=ph_H2O x=logSB;
run;


proc sgplot data=geo;
scatter y=ph_H2O x=mAC;
loess y=ph_H2O x=mAC;
run;

proc tpspline data=geo;
model ph_H2O = (logSB);
output out=a pred;
run;

*proc contents data=a; run;
proc sort data=a; by logSB; run;
proc sgplot data=a;
series y=P_PH_H2O x=logSB;
run;

proc tpspline data=geo;
model ph_H2O = (mAC);
output out=b pred;
run;


*proc contents data=b; run;
proc sort data=b; by mAC; run;
proc sgplot data=b;
series y=P_PH_H2O x=mAC;
run;

/*
data sasuser.geo_MOL_logSB; set a; run;
data sasuser.geo_MOL_mAC; set b; run;
*/

quit;