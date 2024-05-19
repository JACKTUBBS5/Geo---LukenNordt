options nodate nonumber ps=200 ls=80;
title 'SoilGeoChem';


LIBNAME ldata '/home/jacktubbs/my_shared_file_links/jacktubbs/myfolders/Gary_Geology';
*FILENAME user '/home/jacktubbs/my_shared_file_links/jacktubbs/myfolders/Gary_Geology';
/*\*/*\*/*\*/*\*/*\*/*\*/*\*/*\*/*\*/*\*/*\*/*\*/*\*/*\*/*\*/*\*/*\*/*\*


*Read in PPM1 data;
DATA geo;
	SET ldata.PPM1_data;
RUN;


FILENAME REFFILE '/home/jacktubbs/my_shared_file_links/jacktubbs/myfolders
/Gary_Geology/User_data.xls';

PROC IMPORT DATAFILE=REFFILE
	DBMS=XLS
	OUT=WORK.IMPORT;
	GETNAMES=YES;
RUN;


data user; set work.import; run;

*Molal conversions;
DATA user;
	SET user;

	Fe2O3_mol = Fe2O3*10/159.69;
	MnO_mol = MnO*10/70.94;
	P2O5_mol = P2O5*10/283.89;
	SiO2_mol = SiO2*10/60.08;
	TiO2_mol = TiO2*10/79.87;
	ZrO2_mol = ZrO2*10/ 123.22;
	Al2O3_mol = Al2O3*10/101.96;
	CaO_mol = CaO*10/56.08;
	Na2O_mol = Na2O*10/61.98;
	MgO_mol = MgO*10/40.30;
	K2O_mol = K2O*10/94.20;
RUN;

proc sort data=user;
	by Pedon_ID; 
run;


*Combine the data sets;
data combined;
	length PEDON_ID $ 25;
	format Fe2O3 MnO P2O5 SiO2 TiO2 ZrO2 Al2O3 CaO Na2O MgO K2O 
			Fe2O3_mol MnO_mol P2O5_mol SiO2_mol TiO2_mol ZrO2_mol 
			Al2O3_mol CaO_mol Na2O_mol MgO_mol K2O_mol lCaO
			8.2;
	set geo user;

	lCaO = log(CaO);
	lCaO_mol = log(CaO_mol);
run;

data combined; set combined;
if Fe2O3_mol = 0 then Fe2O3_mol = .;
if MnO_mol = 0 then MnO_mol = .;
if P2O5_mol = 0 then P2O5_mol = .;
if SiO2_mol = 0 then SiO2_mol = .;
if TiO2_mol = 0 then TiO2_mol = .;
if ZrO2_mol = 0 then ZrO2_mol = .;
if Al2O3_mol = 0 then Al2O3_mol = .;
if CaO_mol = 0 then CaO_mol = .;
if Na2O_mol = 0 then Na2O_mol = .;
if MgO_mol = 0 then MgO_mol = .;
if K2O_mol = 0 then K2O_mol = .;

run;

*Model;																				
%let MolModel = Fe2O3_mol MnO_mol P2O5_mol SiO2_mol TiO2_mol ZrO2_mol Al2O3_mol 
CaO_mol Na2O_mol MgO_mol K2O_mol;


*-------------------------- Universal model -------------------------------*;		
proc pls data=combined nfac=4 details missing=avg;
	model cl_MAP cl_MAT = &MolModel;
	output out=pls PRESS=PRESS XSCORE=factor;
	title2 'PLS';
run;
proc sort data=pls out=pls; by pedon_id; run;
data user2;
	merge user(in=a) pls;
	by pedon_id;
	if a;
run;

/*
data ldata.pls; set pls; run;
*/

proc tpspline data=PLS;
	model cl_MAT cl_MAP = (factor1 factor2 factor3 factor4);
	where cl_MAT ne . and cl_MAP ne .;
	score data=User2 out=user_pred;
	output out=estimated pred uclm lclm r;
run;


*------------- Universal model adjustments --------------------*;					
data pred_adj;
	set user_pred;
	* subtract 2.98 degrees to p_cl_MAT for Subglacial till;
	if SubglacialTill = 1 then do;
		p_cl_MAT = p_cl_MAT - 2.98;
	end;
	* subtract 1.79 degrees to p_cl_MAT for Till;
	if Till = 1 then do;
		p_cl_MAT = p_cl_MAT - 1.79;
	end;
run;

proc sgplot data=estimated;
scatter y=p_cl_MAT x=cl_MAT;
loess y=p_cl_MAT x=cl_MAT;
run;

proc quantreg data=estimated
                   plots=fitplot ci=none algorithm=interior(tolerance=1e-4) ;
   model p_cl_MAT =  cl_MAT 
             / nosummary quantile=0.1,0.5,0.90;
             output out=q_MAT p=pred q=quant res=r;
run;

data q_MAT; set q_MAT; n=_n_;  
q_10 = pred1;q_50=pred2; q_90=pred3; 
keep Pedon_ID cl_MAT p_cl_MAT q_10 q_50 q_90 n; run;
proc print data=q_MAT; where n <11; run;

proc sgplot data=estimated;
scatter y=p_cl_MAP x=cl_MAP;
loess y=p_cl_MAP x=cl_MAP;
run;

proc quantreg data=estimated
                   plots=fitplot ci=none algorithm=interior(tolerance=1e-4) ;
   model p_cl_MAP =  cl_MAP 
                 / nosummary quantile=0.1,0.5,0.90;
             output out=q_MAT p=pred q=quant res=r;                
run;

data q_MAP; set q_MAP; n=_n_;  
q_10 = pred1;q_50=pred2; q_90=pred3; 
keep Pedon_ID cl_MAP p_cl_MAP q_10 q_50 q_90 n; run;
proc print data=q_MAP; where n <11; run;


*-------------------- MACRO to compute intervals ----------------------;
%macro bounds;
options nonotes;

*Calculate number of user observations;
proc sql noprint;
	select count(*)
	into :nobs
	from Work.user;
quit;

*Create SS (Sums of Squares) variable;
data estimated;
	set estimated;
	SS = factor1**2 + factor2**2 + factor3**2 + factor4**2;
run;

*Merge user input and predicted values;
data user2;
	merge user user_pred;
run;

%do i = 1 %to &nobs;

	*Select the i-th observation from the user data;
	data temp;
		set user2(firstobs = &i);
	run;
	data temp;
		set temp(obs=1);
		SSu = factor1**2 + factor2**2 + factor3**2 + factor4**2;
	run;

	*Store SSu globally as SSobs;
	proc sql noprint;
		select SSu
			into :SSobs
		from Work.temp;
	quit;

	*Merge with data set;
	data compare;
		merge geo estimated;
		SSu = &SSobs;
		diff = abs(SS - SSu);
	run;

	*Sort in order (first obs will be closest spatially);
	proc sort data=compare out=compare; by diff; run;
	
	*Keep the closest match;
	data bounds;
		set compare(obs = 1);
		l_mat = p_cl_mat - lclm_cl_mat;		*distance below mat;
		u_mat = uclm_cl_mat - p_cl_mat;		*distance above mat;
		l_map = p_cl_map - lclm_cl_map;		*distance below mat;
		u_map = uclm_cl_map - p_cl_map;		*distance above mat;
	run;
	
	*Store bounds globally;
	proc sql noprint;
		select l_mat, u_mat, l_map, u_map
			into :l_mat, :u_mat, :l_map, :u_map
		from Work.bounds;
	quit;

	*Calculate interval estimate;
	data interval;
		set user_pred(firstobs = &i);
	run;
	data interval;
		set interval(obs = 1);
		low_MAP = max(0,round(p_cl_map - &l_map, 1));
		high_MAP = round(p_cl_map + &u_map, 1);
		low_MAT = round(max(p_cl_mat - &l_mat, 0), 0.1);
		high_MAT = round(p_cl_mat + &u_mat, 0.1);
	run;

	*Combine iteratively;
	proc append base=predictions data=interval; run;

	%end;

*Rename predicted values;
data predictions;
	merge predictions user(keep=Pedon_ID);
	best_MAP = round(p_cl_MAP,1);
	best_MAT = round(p_cl_MAT,0.1);
	*drop p_cl_MAP p_cl_MAT;
run;

%mend;

*Run the macro;
%bounds;


*Output geology predictions;
data geo_pred;
	merge geo estimated;
run;


*-------------- Print results -------------------*;
title2 'Predictions';
proc print data=predictions;
	var Pedon_ID low_MAP best_MAP high_MAP low_MAT best_MAT high_MAT;
run;

/*
*Kill the library;
proc datasets library=work kill nolist; run; 
quit;
*/






