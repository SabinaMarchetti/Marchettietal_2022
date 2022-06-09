#include <stdio.h>
#include <stdlib.h>
#include "def_.h"
#include <math.h>

int main(int argc, char *argv[])
{
    int t;
	int j=0;
    double *T;
	int i;
	int tt;
	int a,aa;
	int r;
	int radj;
	int* prog;
	double popT;
    double ***Y;
    double *y;
	double h=0.115;
    double* I_1;
	double* I_2;
	double* I_br;
    double* lam_wild;
    double* lam_en;
    double* lam_br;
    double* lambda_wild;
	double* lambda_wildS;
    double* lambda_en;
    double* lambda_br;
	double prev_variant1;
	double II;
	double propFirstDose;
	double policy;
	double HospBTI;
    Param param;
    char polfile[1000];
    char paramfile[1000];
    char popRfile[1000];
    char vaccfile[1000];
    char inffile[1000];
    char seropfile[1000];
    char prevvariant1file[1000];
	char regefffile[1000];
	char mobilfile[1000];
	char sevfile[1000];
	char ifrfile[1000];
	char hfile[1000];
	char ICUfile[1000];
	char proppopfile[1000];
	char contactfile[1000];
	char schoolfile[1000];
	char cov1file[1000];
	char cov2file[1000];
	char cov3file[1000];
	char vcomp1file[1000];
	char vcomp2file[1000];
	char Hpeoplefile[1000];
	char ICUpeoplefile[1000];
	char RemBfile[1000];
	char suscAgefile[1000];
	char temperature1file[1000];
	char temperature2file[1000];
	char BTInffile[1000];
	char efficacy1file[1000];
	char efficacy2file[1000];
	char FirstDosefile[1000];
	char propBTIfile[1000];
	
	/*
	AGE CLASSES:
	0-12
	13-18
	19-64
	65-79
	80+
	*/
	
    sprintf(paramfile,"input/param");
    read_param(paramfile, &param);
    T=(double*) calloc(DAYS, sizeof(double));

	// INPUT FILES
	// ------------------------------------------
	// Static parametrization
	sprintf(proppopfile,"input/proppop");
	sprintf(contactfile,"input/contact");
	sprintf(regefffile,"input/regeff");
	sprintf(suscAgefile,"input/suscAge");
	sprintf(mobilfile,"input/mobil");
	sprintf(sevfile,"input/sev");
	sprintf(hfile,"input/hosp");
	sprintf(ifrfile,"input/ifr");
	sprintf(ICUfile,"input/ICU");
    sprintf(vaccfile,"input/vacc");

	// Dynamic parametrization
	sprintf(polfile, "%s/pol", argv[1]);
	sprintf(cov1file,"%s/cov1",argv[1]);
	sprintf(cov2file,"%s/cov2",argv[1]);
	sprintf(cov3file,"%s/cov3",argv[1]);
	sprintf(schoolfile,"%s/School", argv[1]);
	sprintf(temperature1file,"%s/Temperatures1",argv[1]);
	sprintf(temperature2file,"%s/Temperatures2",argv[1]);
	
	// Model initialization
	sprintf(popRfile,"%s/popR",argv[1]);	    
	sprintf(inffile,"%s/inf", argv[1]);
    sprintf(seropfile,"%s/serop", argv[1]);
    sprintf(prevvariant1file,"%s/prevvariant1", argv[1]);
	sprintf(vcomp1file,"%s/vcomp1",argv[1]);
	sprintf(vcomp2file,"%s/vcomp2",argv[1]);
	sprintf(ICUpeoplefile,"%s/ICUpeople",argv[1]);
	sprintf(Hpeoplefile,"%s/Hpeople",argv[1]);
	sprintf(RemBfile,"%s/RemB",argv[1]);
	sprintf(BTInffile,"%s/BTInf",argv[1]);	
	sprintf(efficacy1file,"%s/efficacy1",argv[1]);
	sprintf(efficacy2file,"%s/efficacy2",argv[1]);
	sprintf(FirstDosefile,"%s/FirstDose",argv[1]);
	sprintf(propBTIfile,"%s/propBTI",argv[1]);
	// ------------------------------------------

	// READ ALL
	// ------------------------------------------
    read_pol(polfile, &param);
	read_vacc(vaccfile, &param);
	read_proppop(proppopfile, &param);
	read_suscAge(suscAgefile, &param);
	read_popreg(popRfile, &param);
    read_inf(inffile, &param);
	read_serop(seropfile, &param);
    read_prevvariant1(prevvariant1file, &param);
	read_regeff(regefffile, &param);
	read_mobil(mobilfile, &param);
	read_sev(sevfile, &param);
	read_ifr(ifrfile, &param);
	read_h(hfile, &param);
	read_ICU(ICUfile, &param);
	read_Hpeople(Hpeoplefile, &param);
	read_ICUpeople(ICUpeoplefile, &param);
	read_contact(contactfile, &param);
	read_school(schoolfile, &param);
	read_coverage_dose1(cov1file, &param);
	read_coverage_dose2(cov2file, &param);
	read_coverage_boost(cov3file, &param);
	read_vcomp1(vcomp1file, &param);
	read_vcomp2(vcomp2file, &param);	
	read_RemB(RemBfile, &param);	
	read_temp_1(temperature1file, &param);
	read_temp_2(temperature2file, &param);
	read_BTInf(BTInffile, &param);
	
	read_efficacy1(efficacy1file, &param);
	read_efficacy2(efficacy2file, &param);
	read_propFirstDose(FirstDosefile, &param);
	read_propBTIage(propBTIfile, &param);
	// ------------------------------------------

    n_eq=17;
	y=(double*)calloc(n_eq,sizeof(double));
	
	prog=(int *)calloc(REGIONS,sizeof(double));
	
	I_1=(double *)calloc(AGECL,sizeof(double));
	I_2=(double *)calloc(AGECL,sizeof(double));
	I_br=(double *)calloc(AGECL,sizeof(double));

    Y=(double ***) calloc(REGIONS, sizeof(double*));	
	
	for(r=0;r<REGIONS;r++){
		Y[r]=(double **) calloc(AGECL, sizeof(double*));
		
		for(a=0;a<AGECL;a++){
			Y[r][a]=(double *) calloc(n_eq, sizeof(double));
			prev_variant1 = param.prev_variant1[r];
			Y[r][a][0]=param.pop[r]*param.proppop[r][a]-(param.inf[r][a]+param.Hpeople[r][a]+param.ICUpeople[r][a]+param.R[r][a]+param.vcomp1[r][a]+param.vcomp2[r][a]+param.RemB[r][a]+param.BTInf[r][a]);
			// Variant 1
			Y[r][a][1]=param.inf[r][a]*prev_variant1*(1.0-0.1*max(param.sev[a],.99)*param.h[a]);
			// Variant 2
			Y[r][a][2]=param.inf[r][a]*(1-prev_variant1)*(1.0-0.1*max(param.sev[a]*param.v2_sev,.99)*param.h[a]);
			// Sympt
			Y[r][a][3]=param.inf[r][a]*0.1*max(param.sev[a]*param.v2_sev,.99)*param.h[a];
			// Hosp
			Y[r][a][4] = param.Hpeople[r][a]*(1-param.propBTIage[a]);
			// ICU
			Y[r][a][5] = param.ICUpeople[r][a]*(1-param.propBTIage[a]);
			//Rec W+E
			Y[r][a][6]=param.R[r][a];
			// First dose
			Y[r][a][7]=param.vcomp1[r][a]*param.propFirstDose[a];
			Y[r][a][8]=param.vcomp2[r][a]*param.propFirstDose[a];
			// Second dose
			Y[r][a][9]=param.vcomp1[r][a]*(1-param.propFirstDose[a]);
			Y[r][a][10]=param.vcomp2[r][a]*(1-param.propFirstDose[a]);
			// BTI wild
			Y[r][a][11]=param.BTInf[r][a]*prev_variant1*max(1-0.4*(1-param.filter)*max(param.sev[a]*param.v2_sev,.99)*param.h[a],1-DUMMY);
			// BTI variant 1
			Y[r][a][12]=param.BTInf[r][a]*(1-prev_variant1)*max(1-0.4*(1-param.filter)*param.sev[a]*param.h[a],1-DUMMY);
			// Symp staging BTI
			Y[r][a][13]=param.BTInf[r][a]*min(0.6*param.filter*max(param.sev[a]*param.v2_sev,.99)*param.h[a],DUMMY);
			// BTI hosp
			Y[r][a][14]= param.Hpeople[r][a]*param.propBTIage[a];
			// BTI ICU
			Y[r][a][15]=param.ICUpeople[r][a]*param.propBTIage[a];
			// BTREM
			Y[r][a][16]=param.RemB[r][a];
			}
	}


    for(t=0; t<DAYS;t++){		
	for(r=0; r<REGIONS; r++){
		param.r=r;
		if((r+1)<4){
			radj = r+1;
			}
		else{
			radj = r+2;
			}
		param.pop[r]=0;
		for(a=0;a<AGECL;a++){
			param.a = a;	
			I_1[a]=Y[r][a][1]+param.contrib*Y[r][a][11];
			
			I_2[a]=Y[r][a][2]+param.contrib*Y[r][a][12];

			for(i=0; i<n_eq; i++){
				param.pop[r]+=Y[r][a][i];
				}
			}
			
		if(t>0){
			tt = 1+(t+1)/7;
			}
		else{
			tt = 1;
			}

		if((((t+1)%7)==0) & (param.pol[r][tt]!=param.pol[r][tt-1])){
			prog[r]=1;
			}
		else{
			if((t==0)){
				if(param.pol[r][tt]==param.pol[r][tt-1]){
				prog[r]=8;
				}
				else{
					prog[r]=1;
				}
			}
			if(prog[r]<10){
				prog[r]+=1;
			}
		}
		
		for(aa=0;aa<AGECL;aa++){
		compute_lambda_1(param.pol[r][tt],param.pol[r][tt-1],prog[r],I_1,tt-1, &param);
		compute_lambda_2(param.pol[r][tt],param.pol[r][tt-1],prog[r],I_2, tt-1,&param);
			}


		for(a=0;a<AGECL;a++){
			for(i=0; i<n_eq; i++){
				// update initial conditions
				y[i]=Y[r][a][i];
				}

			param.a = a;
			param.r = r;
			*T=t;
			simulate(T, t, &param, y);
			II=param.lambda_1[a]*Y[r][a][0]+param.lambda_1[a]*((1-param.v1_effV1)*Y[r][a][7]+(1-param.v2_effV1)*Y[r][a][8]+(1-param.efficacy1[r][a])*Y[r][a][9]+(1-param.efficacy2[r][a])*Y[r][a][10]);
			if(II<0){
				II=0;
			}
			printf("%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t", radj,t, a,param.pop[r]*param.proppop[r][a],param.ifr[a]*(param.hospTime*(Y[r][a][4]+Y[r][a][14])+param.ICUtime*(Y[r][a][5]+Y[r][a][15])),II+param.lambda_2[a]*(Y[r][a][0]+(1-param.v1_effV2)*Y[r][a][7]+(1-param.v2_effV2)*Y[r][a][8]+(1-param.efficacy1[r][a])*Y[r][a][9]+(1-param.efficacy2[r][a])*Y[r][a][10]),II,0.33*(Y[r][a][3]+Y[r][a][13]));
			for(i=0; i<n_eq; i++){
				printf("%lf\t", y[i]);
				}
			printf("\n");

			for(i=0;i<n_eq;i++){
				Y[r][a][i]=y[i];
			}
		}
		}

	}
	exit(0);
	}


