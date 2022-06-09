#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "def_.h"

extern gsl_rng *RNG;

// ----------------------------------------------------------------
// READ PARAMETERS FROM DYNAMIC AND STATIC INPUT FOLDERS
void read_param(char *file, Param *param){
    FILE *fp;
    char str[99];
    int s=0;
    
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    while (fscanf(fp, "%s", str)!=EOF){
        switch(s){
            case 0:
                param->beta_wild=atof(str);
                s++;
                break;
            case 1:
                param->gamma=atof(str);
                s++;
                break;
            case 2:
                param->impr_1=atof(str);
                s++;
                break;
            case 3:
                param->impr_2=atof(str);
                s++;
                break;
            case 4:
                param->v2_sev=atof(str);
                s++;
                break;
            case 5:
                param->ICUtime=atof(str);
                s++;
                break;
            case 6:
                param->hospTime=atof(str);
                s++;
                break;
            default:
                fprintf(stderr, "Impossible param!\n");
                exit(1);
        }
    }
    fclose(fp);
}

void read_pol(char *file, Param *param){
    FILE *fp;
    char str[99];
    int i,j,s;
	int WEEKS;
	
	WEEKS = 1+DAYS/7;

    param->pol=(int **) calloc(REGIONS, sizeof(double *));

    for(i=0; i<REGIONS; i++){
        param->pol[i]=(int *) calloc(WEEKS, sizeof(double));
    }
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file pol for reading: %s\n", strerror(errno));
        exit(1);
    }   
    
    i=0;
    j=0;
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<REGIONS)){
		param->pol[i][j]=atoi(str);
        s++;
        i=s/WEEKS;
        j=s%WEEKS;
    }
    fclose(fp);
}

void read_vacc(char *file, Param *param){
    FILE *fp;
    char str[99];
    int s=0;
    
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    while (fscanf(fp, "%s", str)!=EOF){
        switch(s){
            case 0:
                param->v1_effV1=atof(str);
                s++;
                break;
            case 1:
                param->v2_effV1=atof(str);
                s++;
                break;
            case 2:
                param->v1_effV2=atof(str);
                s++;
                break;
            case 3:
                param->v2_effV2=atof(str);
                s++;
                break;
			case 4:
				param->filter=atof(str);
				s++;
				break;
			case 5:
				param->contrib=atof(str);
				s++;
				break;				
				
            default:
                fprintf(stderr, "Impossible vacc!\n");
                exit(1);
        }
    }
    fclose(fp);
}

void read_proppop(char *file, Param *param){
    FILE *fp;
    char str[99];
    //double str;//[99];
	int s=0,i,j;
	// regionale
    param->proppop=(double **) calloc(REGIONS, sizeof(double));

	for(i=0;i<REGIONS;i++){
		param->proppop[i] = (double *) calloc(AGECL,sizeof(double));
	}
	
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
	i=0;
	j=0;
	s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<REGIONS)){
		param->proppop[i][j]=atof(str);
        s++;
        i=s/AGECL;
        j=s%AGECL;
    }
	
	fclose(fp);	
}

void read_suscAge(char *file, Param *param){
    FILE *fp;
    char str[99];
	int s=0,i=0,j=0;
    param->suscAge=(double **) calloc(AGECL, sizeof(double));

	for(i=0; i<AGECL; i++){
        param->suscAge[i]=(double *) calloc(2, sizeof(double));
    }

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    i=0;
    j=0;
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (j<2)){
		param->suscAge[i][j]=atof(str);
        s++;
        i=s/2;
        j=s%2;
    }
    fclose(fp);
}

void read_popreg(char *file, Param *param){
    FILE *fp;
    char str[99];
	int s;

    param->pop=(double *) calloc(REGIONS, sizeof(double));
	
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }

    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (s<REGIONS)){
		param->pop[s]=atof(str);
		s++;
		}

    fclose(fp);	
}

void read_inf(char *file, Param *param){
    FILE *fp;
    char str[99];
    //double str;//[99];
	int i,j,s;
	// regionale
    param->inf=(double **) calloc(REGIONS, sizeof(double*));
	
	for(i=0;i<REGIONS;i++){
		param->inf[i] = (double *) calloc(AGECL,sizeof(double));
	}
	
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
	i=0;
	j=0;
	s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<REGIONS)){
		param->inf[i][j]=atof(str);
        s++;
        i=s/AGECL;
        j=s%AGECL;
    }
    fclose(fp);	
}

void read_serop(char *file, Param *param){
    FILE *fp;
    char str[99];
	int i,j,s;
    param->R=(double **) calloc(REGIONS, sizeof(double*));
	
	for(i=0;i<REGIONS;i++){
		param->R[i] = (double *) calloc(AGECL,sizeof(double));
	}
	
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
	i=0;
	j=0;
	s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<REGIONS)){
		param->R[i][j]=atof(str);
        s++;
        i=s/AGECL;
        j=s%AGECL;
    }
    fclose(fp);
}

void read_prevvariant1(char *file, Param *param){
    FILE *fp;
    char str[99];
	int s;
    param->prev_variant1=(double *) calloc(REGIONS, sizeof(double));

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (s<REGIONS)){
		param->prev_variant1[s]=atof(str);
        s++;
    }

    fclose(fp);
}

void read_regeff(char *file, Param *param){
    FILE *fp;
    char str[99];
    //double str;//[99];
	int s=0;
	// regionale
    param->regeff=(double *) calloc(REGIONS, sizeof(double));

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    while ((fscanf(fp, "%s", str)!=EOF) & (s < REGIONS)){
        param->regeff[s]=atof(str);		
        s++;
    }
    fclose(fp);
}

void read_mobil(char *file, Param *param){
    FILE *fp;
    char str[99];
	int i,j;
	int s=0;

    param->mobil=(double **) calloc(COLORS+1, sizeof(double *));

    for(i=0; i<(COLORS+1); i++){
        param->mobil[i]=(double *) calloc(2, sizeof(double));
    }
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    i=0;
	j=0;
	s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i < (COLORS+1))){
        param->mobil[i][j]=atof(str);
        s++;
		i=s/2;
		j=s%2;
    }

    fclose(fp);
	}

void read_sev(char *file, Param *param){
    FILE *fp;
    char str[99];
	int s;
    param->sev=(double *) calloc(AGECL, sizeof(double));

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    while ((fscanf(fp, "%s", str)!=EOF) & (s < AGECL)){
        param->sev[s]=atof(str);			
        s++;
    }
    fclose(fp);
}

void read_ifr(char *file, Param *param){
    FILE *fp;
    char str[99];
    //double str;//[99];
	int s=0;
	// regionale
    param->ifr=(double *) calloc(AGECL, sizeof(double));

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    while ((fscanf(fp, "%s", str)!=EOF) & (s < AGECL)){
        param->ifr[s]=atof(str);
        s++;
    }
    fclose(fp);
}

void read_h(char *file, Param *param){
    FILE *fp;
    char str[99];
	int s=0;

    param->h=(double *) calloc(AGECL, sizeof(double));

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    while ((fscanf(fp, "%s", str)!=EOF) & (s < AGECL)){
        param->h[s]=atof(str);
        s++;
    }
    fclose(fp);
}

void read_ICU(char *file, Param *param){
    FILE *fp;
    char str[99];
    //double str;//[99];
	int s=0;
	// regionale
    param->ICU=(double *) calloc(AGECL, sizeof(double));

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    while ((fscanf(fp, "%s", str)!=EOF) & (s < AGECL)){
        param->ICU[s]=atof(str);		
        s++;
    }
    fclose(fp);
}

void read_Hpeople(char *file, Param *param){
    FILE *fp;
    char str[99];
    //double str;//[99];
	int i,j,s;
	// regionale
    param->Hpeople=(double **) calloc(REGIONS, sizeof(double*));
	
	for(i=0;i<REGIONS;i++){
		param->Hpeople[i] = (double *) calloc(AGECL,sizeof(double));
	}
	
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
	i=0;
	j=0;
	s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<REGIONS)){
		param->Hpeople[i][j]=atof(str);
        s++;
        i=s/AGECL;
        j=s%AGECL;
    }
    fclose(fp);	
}

void read_ICUpeople(char *file, Param *param){
    FILE *fp;
    char str[99];
    //double str;//[99];
	int i,j,s;
	// regionale
    param->ICUpeople=(double **) calloc(REGIONS, sizeof(double*));
	
	for(i=0;i<REGIONS;i++){
		param->ICUpeople[i] = (double *) calloc(AGECL,sizeof(double));
	}
	
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
	i=0;
	j=0;
	s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<REGIONS)){
		param->ICUpeople[i][j]=atof(str);
        s++;
        i=s/AGECL;
        j=s%AGECL;
    }
    fclose(fp);	
}

void read_contact(char *file, Param *param){
    FILE *fp;
    char str[99];
	int i,j,s;

    param->C=(double **) calloc(AGECL, sizeof(double*));
	
	for(i=0; i<REGIONS; i++){
        param->C[i]=(double *) calloc(AGECL, sizeof(double));
    }

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    i=0;
    j=0;
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<AGECL)){
		param->C[i][j]=atof(str);
        s++;
        i=s/AGECL;
        j=s%AGECL;
    }
    fclose(fp);
}

void read_school(char *file, Param *param){
    FILE *fp;
    char str[99];
	int i,j,s;
	
    param->School=(double **) calloc(REGIONS, sizeof(double *));

    for(i=0; i<REGIONS; i++){
        param->School[i]=(double *) calloc(3, sizeof(double));
    }
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file for reading: %s\n", strerror(errno));
        exit(1);
    }   
    
    i=0;
    j=0;
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<REGIONS)){
		param->School[i][j]=atof(str);
        s++;
        i=s/3;
        j=s%3;
    }
    fclose(fp);
}

void read_coverage_dose1(char *file, Param *param){
    FILE *fp;
    char str[99];
	int i,j,s;
	
    param->cov1=(double **) calloc(AGECL, sizeof(double*));
	
	for(i=0; i<REGIONS; i++){
        param->cov1[i]=(double *) calloc(REGIONS, sizeof(double));
    }

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    i=0;
    j=0;
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<AGECL)){
		param->cov1[i][j]=atof(str)/30;
        s++;
        i=s/REGIONS;
        j=s%REGIONS;
    }
    fclose(fp);
}

void read_coverage_dose2(char *file, Param *param){
    FILE *fp;
    char str[99];
	int i,j,s;
    param->cov2=(double **) calloc(AGECL, sizeof(double*));
	
	for(i=0; i<REGIONS; i++){
        param->cov2[i]=(double *) calloc(REGIONS, sizeof(double));
    }

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    i=0;
    j=0;
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<AGECL)){
		param->cov2[i][j]=atof(str)/30;
		s++;
        i=s/REGIONS;
        j=s%REGIONS;
    }
    fclose(fp);
}

void read_coverage_boost(char *file, Param *param){
    FILE *fp;
    char str[99];
    //double str;//[99];
	int i,j,s;
	// regionale
    param->cov3=(double **) calloc(AGECL, sizeof(double*));
	
	for(i=0; i<REGIONS; i++){
        param->cov3[i]=(double *) calloc(REGIONS, sizeof(double));
    }

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    i=0;
    j=0;
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<AGECL)){
		param->cov3[i][j]=atof(str);
        s++;
        i=s/REGIONS;
        j=s%REGIONS;
    }
    fclose(fp);
}

void read_vcomp1(char *file, Param *param){
    FILE *fp;
    char str[99];
    //double str;//[99];
	int i,j,s;
	// regionale
    param->vcomp1=(double **) calloc(REGIONS, sizeof(double*));
	
	for(i=0;i<REGIONS;i++){
		param->vcomp1[i] = (double *) calloc(AGECL,sizeof(double));
	}
	
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
	i=0;
	j=0;
	s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<REGIONS)){
		param->vcomp1[i][j]=atof(str);
        s++;
        i=s/AGECL;
        j=s%AGECL;
    }
    fclose(fp);	
}

void read_vcomp2(char *file, Param *param){
    FILE *fp;
    char str[99];
	int i,j,s;

    param->vcomp2=(double **) calloc(REGIONS, sizeof(double*));
	
	for(i=0;i<REGIONS;i++){
		param->vcomp2[i] = (double *) calloc(AGECL,sizeof(double));
	}
	
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
	i=0;
	j=0;
	s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<REGIONS)){
		param->vcomp2[i][j]=atof(str);
        s++;
        i=s/AGECL;
        j=s%AGECL;
    }
    fclose(fp);	
}

void read_RemB(char *file, Param *param){
    FILE *fp;
    char str[99];
	int i,j,s;
	// regionale
    param->RemB=(double **) calloc(REGIONS, sizeof(double*));
	
	for(i=0;i<REGIONS;i++){
		param->RemB[i] = (double *) calloc(AGECL,sizeof(double));
	}
	
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
	i=0;
	j=0;
	s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<REGIONS)){
		param->RemB[i][j]=atof(str);
        s++;
        i=s/AGECL;
        j=s%AGECL;
    }
    fclose(fp);	
}

void read_temp_1(char *file, Param *param){
    FILE *fp;
    char str[99];
    int i,j,s;
	int WEEKS;
	
	WEEKS = 3;

    param->Temp_1=(double **) calloc(REGIONS, sizeof(double *));

    for(i=0; i<REGIONS; i++){
        param->Temp_1[i]=(double *) calloc(WEEKS, sizeof(double));
    }
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file temp1 for reading: %s\n", strerror(errno));
        exit(1);
    }   
    
    i=0;
    j=0;
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF)&(i<REGIONS)){
		param->Temp_1[i][j]=atof(str);
        s++;
        i=s/WEEKS;
        j=s%WEEKS;
    }

    fclose(fp);
}

void read_temp_2(char *file, Param *param){
    FILE *fp;
    char str[99];
    int i,j,s;
	int WEEKS;
	
	WEEKS = 3;

    param->Temp_2=(double **) calloc(REGIONS, sizeof(double *));

    for(i=0; i<REGIONS; i++){
        param->Temp_2[i]=(double *) calloc(WEEKS, sizeof(double));
    }
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file temp2 for reading: %s\n", strerror(errno));
        exit(1);
    }   
    
    i=0;
    j=0;
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF)&(i<REGIONS)){
		param->Temp_2[i][j]=atof(str);
        s++;
        i=s/WEEKS;
        j=s%WEEKS;
    }

    fclose(fp);
}

void read_BTInf(char *file, Param *param){
    FILE *fp;
    char str[99];
    //double str;//[99];
	int i,j,s;
	// regionale
    param->BTInf=(double **) calloc(REGIONS, sizeof(double*));
	
	for(i=0;i<REGIONS;i++){
		param->BTInf[i] = (double *) calloc(AGECL,sizeof(double));
	}
	
    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
	i=0;
	j=0;
	s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<REGIONS)){
		param->BTInf[i][j]=atof(str);
        s++;
        i=s/AGECL;
        j=s%AGECL;
    }
    fclose(fp);	
}

void read_efficacy1(char *file, Param *param){
    FILE *fp;
    char str[99];
    //double str;//[99];
	int i,j,s;
	
	// 12 x (agecl-Pfizer ; agecl-AZ)
    param->efficacy1=(double **) calloc(REGIONS, sizeof(double*));
	
	for(i=0; i<REGIONS; i++){
        param->efficacy1[i]=(double *) calloc(AGECL, sizeof(double));
    }

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    i=0;
    j=0;
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (j<REGIONS)){
		param->efficacy1[i][j]=atof(str)/100;
        s++;
        i=s/AGECL;
        j=s%AGECL;
    }
	//exit(1);
    fclose(fp);
}

void read_efficacy2(char *file, Param *param){
    FILE *fp;
    char str[99];
    //double str;//[99];
	int i,j,s;
	
	// 12 x (agecl-Pfizer ; agecl-AZ)
    param->efficacy2=(double **) calloc(REGIONS, sizeof(double*));
	
	for(i=0; i<REGIONS; i++){
        param->efficacy2[i]=(double *) calloc(AGECL, sizeof(double));
    }

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    i=0;
    j=0;
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (j<REGIONS)){
		param->efficacy2[i][j]=atof(str)/100;
        s++;
        i=s/AGECL;
        j=s%AGECL;
    }
	//exit(1);
    fclose(fp);
}

void read_propFirstDose(char *file, Param *param){
    FILE *fp;
    char str[99];
	int i,j,s;
	
    param->propFirstDose=(double *) calloc(AGECL, sizeof(double));

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<AGECL)){
		param->propFirstDose[s]=atof(str);
        s++;
    }
    fclose(fp);
}

void read_propBTIage(char *file, Param *param){
    FILE *fp;
    char str[99];
	int i,j,s;
    param->propBTIage=(double *) calloc(AGECL, sizeof(double));

    if((fp=fopen(file,"r"))==NULL){
        fprintf(stderr,"Error opening file %s for reading\n",file);
        exit(1);
    }
    
    s=0;
    while ((fscanf(fp, "%s", str)!=EOF) & (i<AGECL)){
		param->propBTIage[s]=atof(str);
        s++;
    }
    fclose(fp);
}

//-------------------------------------------------------------------------
// COMPUTE FOI

void compute_lambda_1(int restr,int previous,int prog, double* I_1, int t,Param *param){
    int r,i,k;
	double pprog;
	double p;

	
	param->lambda_1=(double *) calloc(AGECL, sizeof(double));

	p=prog*0.1;
	pprog = param->mobil[previous][0]*(1-p)+param->mobil[restr][0]*p;

	r = param->r;
		for(i=0; i<AGECL; i++){
		for (k=0;k<AGECL;k++){
			param->lambda_1[i]+=param->C[i][k]*I_1[k]/(param->pop[r]*param->proppop[r][k]);
			if(k<2){
				param->lambda_1[i]*=param->School[r][t];
			}
			if(k>1 & k<4){
				param->lambda_1[i]*=param->mobil[COLORS][restr];
			}
			}
    param->lambda_1[i]*=param->impr_1*param->suscAge[i][0]*pprog*param->regeff[r]*param->beta_wild*(1+param->Temp_1[r][t]);
		}
		
    }

	

void compute_lambda_2(int restr,int previous,int prog, double* I_2,int t, Param *param){
    int r,i,k;
	double pprog;
	double p;
	
	param->lambda_2=(double *) calloc(AGECL, sizeof(double));

	p=prog*0.1;
	pprog = param->mobil[previous][0]*(1-p)+param->mobil[restr][0]*p;
	r = param->r;
		for(i=0; i<AGECL; i++){
		for (k=0;k<AGECL;k++){
			param->lambda_2[i]+=param->C[i][k]*I_2[k]/(param->pop[r]*param->proppop[r][k]);

			if(k<2){
				param->lambda_2[i]*=param->School[r][t];
			}
			if(k>1 & k<4){
				param->lambda_2[i]*=param->mobil[COLORS][restr];
			}
			}

    param->lambda_2[i]*=param->impr_2*param->suscAge[i][1]*pprog*param->regeff[r]*param->beta_wild*(1+param->Temp_2[r][t]);
        }
    }
	

//--------------------------------------------------------------------------------------------------
// MODEL EQUATIONS

int eq_sys (double t,const double y[], double f[],void *params)
{
    Param *param;
    param=(Param *)params;
	int a=param->a;
	int r=param->r;
	double sigma;
	double boostV1,boostV2;
	double sev;
	double gammavacc;
	
	gammavacc = param->gamma*1.17;
	
	sev=min(param->sev[a]*param->v2_sev,.99);
	
	sigma=0.59;//0.63;
	boostV1=0.02857143; // 35 days**(-1)
	boostV2 =0.01111111; //90 days**(-1)

	// NATURAL DISEASE
	
	// Susceptibles
    f[0] = -(param->lambda_1[a]+param->lambda_2[a]+param->cov1[a][r]+param->cov2[a][r])*y[0];
	// Variant 1
	f[1] = param->lambda_1[a]*y[0]-param->gamma*y[1];
	// Variant 2
	f[2] = param->lambda_2[a]*y[0]-param->gamma*y[2];
	// Delay
	f[3] = param->gamma*param->h[a]*(sev*y[1]+param->sev[a]*y[2])-sigma*y[3];
	// MA
	f[4] = sigma*(1-param->ICU[a])*y[3]-(param->hospTime*1)*y[4];
	// ICU
	f[5]=sigma*param->ICU[a]*y[3]-param->ICUtime*y[5];
	// Recovered 1+2
	f[6] = param->gamma*((1-param->h[a]*param->sev[a])*y[2]+(1-param->h[a]*sev)*y[1])+(1-param->ifr[a])*(param->hospTime*1*y[4]+param->ICUtime*1*y[5])-(param->cov1[a][r]+param->cov2[a][r])*y[6];
	
	// IMMUNIZED
	// Vaccines: First dose
	// First dose Vaccine group 1
	f[7] = param->cov1[a][r]*y[0]-(param->lambda_1[a]*(1-param->v1_effV1)+param->lambda_2[a]*(1-param->v1_effV2)+boostV1)*y[7];
	// First dose Vaccine group 2
	f[8] = param->cov2[a][r]*y[0]-(param->lambda_1[a]*(1-param->v2_effV1)+param->lambda_2[a]*(1-param->v2_effV2)+boostV2)*y[8];
	// Vaccines: Second dose
	// Second dose Vaccine group 1
	f[9] = boostV1*y[7]+param->cov3[a][r]*y[10]-(param->lambda_1[a]*(1-param->efficacy1[r][a])+param->lambda_2[a]*(1-param->efficacy1[r][a]))*y[9];
	// Second dose Vaccine group 2
	f[10] = boostV2*y[8]-(param->lambda_1[a]*(1-param->efficacy2[r][a])+param->lambda_2[a]*(1-param->efficacy2[r][a])+param->cov3[a][r])*y[10];	
	// BTI: Variant 1
	f[11] = param->lambda_1[a]*((1-param->v1_effV1)*y[7]+(1-param->v2_effV1)*y[8]+(1-param->efficacy1[r][a])*y[9]+(1-param->efficacy1[r][a])*y[10])-gammavacc*y[11];    
	// BTI: Variant 2
	f[12] = param->lambda_2[a]*((1-param->v1_effV2)*y[7]+(1-param->v2_effV2)*y[8]+(1-param->efficacy2[r][a])*y[9]+(1-param->efficacy2[r][a])*y[10])-gammavacc*y[12];
	// Delay
	f[13] = gammavacc*param->filter*param->h[a]*(sev*y[11]+param->sev[a]*y[12])-sigma*y[13];
	// MA_BT
	f[14] = sigma*(1-param->ICU[a])*y[13]-param->hospTime*y[14];
	// ICU_BT
	f[15]=sigma*param->ICU[a]*y[13]-param->ICUtime*y[15];
	// Recovered & Vaccinated
	f[16] = gammavacc*((1-param->h[a]*param->sev[a]*param->filter)*y[12]+(1-param->h[a]*sev*param->filter)*y[11])+(1-param->ifr[a])*(param->hospTime*y[14]+param->ICUtime*y[15])+(param->cov1[a][r]+param->cov2[a][r])*y[6];
	return GSL_SUCCESS;
}

// --------------------------------------------------------------------
// INTEGRATION
void simulate(double *T, int t, Param *param, double *y){
    double h;
    int j,s;
    
    const gsl_odeiv2_step_type *Typ = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step *ste = gsl_odeiv2_step_alloc(Typ, n_eq);
    gsl_odeiv2_control *con = gsl_odeiv2_control_y_new(1e-4, 1e-4);
    gsl_odeiv2_evolve *evo = gsl_odeiv2_evolve_alloc(n_eq);
    
    gsl_odeiv2_system sys = { eq_sys, NULL, n_eq, param};
    
    h=0.01;
    double tend=t+(1.0);
    
    while(*T<tend){
        s=gsl_odeiv2_evolve_apply(evo, con, ste, &sys, T, tend, &h, y);
        
        for(j=0; j<n_eq; j++)
            if(y[j]<0)
                y[j]=0;
       
        if (s != GSL_SUCCESS){
            fprintf (stderr, "error: driver returned %d\n", s);
            exit(1);
        }
    }
    
    gsl_odeiv2_evolve_free(evo);
    gsl_odeiv2_control_free(con);
    gsl_odeiv2_step_free(ste);
    return;
}
