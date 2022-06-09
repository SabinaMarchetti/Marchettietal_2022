#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define DELTA 365
#define DAYS 21
#define REGIONS 21
#define AGECL 5
#define COLORS 5
#define max(x,y) (x>y?x:y)
#define min(x,y) (x<y?x:y)
#define DUMMY 1

int n_eq;
double *ag;

typedef struct{
	int a;
	int r;
    int** pol;
	double** mobil;
    double* pop;
    double** inf;
    double** R;
	double** vcomp1;
	double** vcomp2;
    double* prev_variant1;
	double* sev;
	double** C;
	double* regeff;
	double* h;
	double** suscAge;
	double filter;
	double* propBTIage;
    double beta_wild;
	double impr_1;
	double impr_2;
    double gamma;
	double hospTime;
	double** School;
	double* ifr;
	double ICUtime;
    double v1_effV1;
    double v1_effV2;
    double v2_effV1;
	double v2_effV2;
	double v2_sev;
	double temp_mob;
	double** efficacy1;
	double** efficacy2;
	double contrib;
	double gammavacc;
	double* propFirstDose;
	double propBTI;
    double** cov1;
	double** cov2;
	double** cov3;
    double* lambda_1;
    double* lambda_2;
	double* ICU;
	double** proppop;
	double** ICUpeople;
	double** Hpeople;
	double** RemB;
	double** Temp_1;
	double** Temp_2;
	double** BTInf;

} Param;

void compute_lambda_1(int restr,int previous,int prog, double* y, int t,Param *param);
void compute_lambda_2(int restr,int previous,int prog, double* y, int t,Param *param);
void read_pol(char *file, Param *param);
void read_popreg(char *file, Param *param);
void read_contact(char *file, Param *param);
void read_vcomp1(char *file, Param *param);
void read_vcomp2(char *file, Param *param);
void read_coverage_dose1(char *file, Param *param);
void read_coverage_dose2(char *file, Param *param);
void read_coverage_boost(char *file, Param *param);
void read_suscAge(char *file, Param *param);
void read_inf(char *file, Param *param);
void read_param(char *file,Param *param);
void read_vacc(char *file,Param *param);
void read_serop(char *file, Param *param);
void read_prevvariant1(char *file, Param *param);
void read_regeff(char *file, Param *param);
void read_mobil(char *file, Param *param);
void read_sev(char *file, Param *param);
void read_ifr(char *file, Param *param);
void read_h(char *file, Param *param);
void read_ICU(char *file, Param *param);
void read_school(char *file, Param *param);
void read_proppop(char *file, Param *param);
void read_ICUpeople(char *file, Param *param);
void read_Hpeople(char *file, Param *param);
void read_RemB(char *file, Param *param);
void read_temp_1(char *file, Param *param);
void read_temp_2(char *file, Param *param);
void read_BTInf(char *file, Param *param);
void read_efficacy1(char *file, Param *param);
void read_efficacy2(char *file, Param *param);
void read_propFirstDose(char *file, Param *param);
void read_propBTIage(char *file, Param *param);
int eq_sys (double t, const double y[], double f[],void *params);
void simulate(double *T,int t, Param *param, double *y);

gsl_rng *RNG;