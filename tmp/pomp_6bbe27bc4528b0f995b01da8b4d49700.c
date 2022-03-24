/* pomp C snippet file: pomp_6bbe27bc4528b0f995b01da8b4d49700 */
/* Time: 2022-03-22 16:16:22.755 -0400 */
/* Salt: 04DA552A65239912995C4995 */

#include <C:/Users/ntreu/Documents/R/win-library/4.1/pomp/include/pomp.h>
#include <R_ext/Rdynload.h>

 


/* C snippet: 'rinit' */
#define rho		(__p[__parindex[0]])
#define tau		(__p[__parindex[1]])
#define beta1		(__p[__parindex[2]])
#define beta2		(__p[__parindex[3]])
#define beta3		(__p[__parindex[4]])
#define beta4		(__p[__parindex[5]])
#define beta5		(__p[__parindex[6]])
#define beta6		(__p[__parindex[7]])
#define betat		(__p[__parindex[8]])
#define gamma		(__p[__parindex[9]])
#define sigma		(__p[__parindex[10]])
#define theta0		(__p[__parindex[11]])
#define alpha		(__p[__parindex[12]])
#define mu		(__p[__parindex[13]])
#define delta		(__p[__parindex[14]])
#define nu		(__p[__parindex[15]])
#define pop_0		(__p[__parindex[16]])
#define sig_sq		(__p[__parindex[17]])
#define S_0		(__p[__parindex[18]])
#define E_0		(__p[__parindex[19]])
#define I_0		(__p[__parindex[20]])
#define A_0		(__p[__parindex[21]])
#define R_0		(__p[__parindex[22]])
#define seas1		(__covars[__covindex[0]])
#define seas2		(__covars[__covindex[1]])
#define seas3		(__covars[__covindex[2]])
#define seas4		(__covars[__covindex[3]])
#define seas5		(__covars[__covindex[4]])
#define seas6		(__covars[__covindex[5]])
#define S		(__x[__stateindex[0]])
#define E		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define A		(__x[__stateindex[3]])
#define R		(__x[__stateindex[4]])
#define incid		(__x[__stateindex[5]])
#define foival		(__x[__stateindex[6]])
#define Str0		(__x[__stateindex[7]])
#define Sout		(__x[__stateindex[8]])
#define Sin		(__x[__stateindex[9]])

void __pomp_rinit (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)
{
 double frac = pop_0 / (S_0 + E_0 + I_0 + A_0 + R_0); 
 S = nearbyint(frac * S_0); 
 E = nearbyint(frac * E_0); 
 I = nearbyint(frac * I_0); 
 A = nearbyint(frac * A_0); 
 R = nearbyint(frac * R_0); 
 incid = 0.0; 
 foival = 0.0; 
 Str0 = 0.0; 
 Sout = 0.0; 
 Sin = 0.0; 
  
}

#undef rho
#undef tau
#undef beta1
#undef beta2
#undef beta3
#undef beta4
#undef beta5
#undef beta6
#undef betat
#undef gamma
#undef sigma
#undef theta0
#undef alpha
#undef mu
#undef delta
#undef nu
#undef pop_0
#undef sig_sq
#undef S_0
#undef E_0
#undef I_0
#undef A_0
#undef R_0
#undef seas1
#undef seas2
#undef seas3
#undef seas4
#undef seas5
#undef seas6
#undef S
#undef E
#undef I
#undef A
#undef R
#undef incid
#undef foival
#undef Str0
#undef Sout
#undef Sin

/* C snippet: 'step.fn' */
#define rho		(__p[__parindex[0]])
#define tau		(__p[__parindex[1]])
#define beta1		(__p[__parindex[2]])
#define beta2		(__p[__parindex[3]])
#define beta3		(__p[__parindex[4]])
#define beta4		(__p[__parindex[5]])
#define beta5		(__p[__parindex[6]])
#define beta6		(__p[__parindex[7]])
#define betat		(__p[__parindex[8]])
#define gamma		(__p[__parindex[9]])
#define sigma		(__p[__parindex[10]])
#define theta0		(__p[__parindex[11]])
#define alpha		(__p[__parindex[12]])
#define mu		(__p[__parindex[13]])
#define delta		(__p[__parindex[14]])
#define nu		(__p[__parindex[15]])
#define pop_0		(__p[__parindex[16]])
#define sig_sq		(__p[__parindex[17]])
#define S_0		(__p[__parindex[18]])
#define E_0		(__p[__parindex[19]])
#define I_0		(__p[__parindex[20]])
#define A_0		(__p[__parindex[21]])
#define R_0		(__p[__parindex[22]])
#define seas1		(__covars[__covindex[0]])
#define seas2		(__covars[__covindex[1]])
#define seas3		(__covars[__covindex[2]])
#define seas4		(__covars[__covindex[3]])
#define seas5		(__covars[__covindex[4]])
#define seas6		(__covars[__covindex[5]])
#define S		(__x[__stateindex[0]])
#define E		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define A		(__x[__stateindex[3]])
#define R		(__x[__stateindex[4]])
#define incid		(__x[__stateindex[5]])
#define foival		(__x[__stateindex[6]])
#define Str0		(__x[__stateindex[7]])
#define Sout		(__x[__stateindex[8]])
#define Sin		(__x[__stateindex[9]])

void __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t, double dt)
{
 double Srate[2]; 
 double Erate[3]; 
 double Irate[2]; 
 double Arate[2]; 
 double Rrate[2]; 
 double Strans[2]; 
 double Etrans[3]; 
 double Itrans[2]; 
 double Atrans[2]; 
 double Rtrans[2]; 
 
 double dgamma; 
 int pop = S + E + I + A + R; 
 int births = rpois(mu*pop*dt); 
 double mybeta = exp(log(beta1*seas1 + beta2*seas2 + beta3*seas3 + beta4*seas4 + beta5*seas5 + beta6*seas6) - (betat)*((t-215)/(430-215))); 
 double foi = pow(I, nu) * mybeta / pop; 
 
 dgamma = rgammawn(sig_sq, dt); 
 foi = foi * dgamma/dt; 
 Srate[0] = foi; 
 Srate[1] = delta; 
 Erate[0] = sigma*(1-theta0); 
 Erate[1] = sigma*theta0; 
 Erate[2] = delta; 
 Irate[0] = gamma; 
 Irate[1] = delta; 
 Arate[0] = gamma; 
 Arate[1] = delta; 
 Rrate[0] = alpha; 
 Rrate[1] = delta; 
 reulermultinom(2,S,&Srate[0],dt,&Strans[0]); 
 reulermultinom(3,E,&Erate[0],dt,&Etrans[0]); 
 reulermultinom(2,I,&Irate[0],dt,&Itrans[0]); 
 reulermultinom(2,A,&Arate[0],dt,&Atrans[0]); 
 reulermultinom(2,R,&Rrate[0],dt,&Rtrans[0]); 
 S += -Strans[0] - Strans[1] + Rtrans[0] + births; 
 E += -Etrans[0] - Etrans[1] - Etrans[2] + Strans[0]; 
 I += -Itrans[0] - Itrans[1] + Etrans[0]; 
 A += -Atrans[0] - Atrans[1] + Etrans[1]; 
 R += -Rtrans[0] - Rtrans[1] + Itrans[0] + Atrans[0]; 
 foival += foi; 
 Str0 += Strans[0]; 
 Sin += Rtrans[0] + births; 
 incid += Etrans[0]; 
 Sout += Strans[0] + Strans[1]; 
  
}

#undef rho
#undef tau
#undef beta1
#undef beta2
#undef beta3
#undef beta4
#undef beta5
#undef beta6
#undef betat
#undef gamma
#undef sigma
#undef theta0
#undef alpha
#undef mu
#undef delta
#undef nu
#undef pop_0
#undef sig_sq
#undef S_0
#undef E_0
#undef I_0
#undef A_0
#undef R_0
#undef seas1
#undef seas2
#undef seas3
#undef seas4
#undef seas5
#undef seas6
#undef S
#undef E
#undef I
#undef A
#undef R
#undef incid
#undef foival
#undef Str0
#undef Sout
#undef Sin

/* C snippet: 'rmeasure' */
#define rho		(__p[__parindex[0]])
#define tau		(__p[__parindex[1]])
#define beta1		(__p[__parindex[2]])
#define beta2		(__p[__parindex[3]])
#define beta3		(__p[__parindex[4]])
#define beta4		(__p[__parindex[5]])
#define beta5		(__p[__parindex[6]])
#define beta6		(__p[__parindex[7]])
#define betat		(__p[__parindex[8]])
#define gamma		(__p[__parindex[9]])
#define sigma		(__p[__parindex[10]])
#define theta0		(__p[__parindex[11]])
#define alpha		(__p[__parindex[12]])
#define mu		(__p[__parindex[13]])
#define delta		(__p[__parindex[14]])
#define nu		(__p[__parindex[15]])
#define pop_0		(__p[__parindex[16]])
#define sig_sq		(__p[__parindex[17]])
#define S_0		(__p[__parindex[18]])
#define E_0		(__p[__parindex[19]])
#define I_0		(__p[__parindex[20]])
#define A_0		(__p[__parindex[21]])
#define R_0		(__p[__parindex[22]])
#define seas1		(__covars[__covindex[0]])
#define seas2		(__covars[__covindex[1]])
#define seas3		(__covars[__covindex[2]])
#define seas4		(__covars[__covindex[3]])
#define seas5		(__covars[__covindex[4]])
#define seas6		(__covars[__covindex[5]])
#define S		(__x[__stateindex[0]])
#define E		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define A		(__x[__stateindex[3]])
#define R		(__x[__stateindex[4]])
#define incid		(__x[__stateindex[5]])
#define foival		(__x[__stateindex[6]])
#define Str0		(__x[__stateindex[7]])
#define Sout		(__x[__stateindex[8]])
#define Sin		(__x[__stateindex[9]])
#define cases		(__y[__obsindex[0]])

void __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 cases = rnbinom_mu(tau, rho*incid);
    if (cases > 0.0) {
      cases = nearbyint(cases);
    } else {
      cases = 0.0;
    }
   
}

#undef rho
#undef tau
#undef beta1
#undef beta2
#undef beta3
#undef beta4
#undef beta5
#undef beta6
#undef betat
#undef gamma
#undef sigma
#undef theta0
#undef alpha
#undef mu
#undef delta
#undef nu
#undef pop_0
#undef sig_sq
#undef S_0
#undef E_0
#undef I_0
#undef A_0
#undef R_0
#undef seas1
#undef seas2
#undef seas3
#undef seas4
#undef seas5
#undef seas6
#undef S
#undef E
#undef I
#undef A
#undef R
#undef incid
#undef foival
#undef Str0
#undef Sout
#undef Sin
#undef cases

/* C snippet: 'dmeasure' */
#define rho		(__p[__parindex[0]])
#define tau		(__p[__parindex[1]])
#define beta1		(__p[__parindex[2]])
#define beta2		(__p[__parindex[3]])
#define beta3		(__p[__parindex[4]])
#define beta4		(__p[__parindex[5]])
#define beta5		(__p[__parindex[6]])
#define beta6		(__p[__parindex[7]])
#define betat		(__p[__parindex[8]])
#define gamma		(__p[__parindex[9]])
#define sigma		(__p[__parindex[10]])
#define theta0		(__p[__parindex[11]])
#define alpha		(__p[__parindex[12]])
#define mu		(__p[__parindex[13]])
#define delta		(__p[__parindex[14]])
#define nu		(__p[__parindex[15]])
#define pop_0		(__p[__parindex[16]])
#define sig_sq		(__p[__parindex[17]])
#define S_0		(__p[__parindex[18]])
#define E_0		(__p[__parindex[19]])
#define I_0		(__p[__parindex[20]])
#define A_0		(__p[__parindex[21]])
#define R_0		(__p[__parindex[22]])
#define seas1		(__covars[__covindex[0]])
#define seas2		(__covars[__covindex[1]])
#define seas3		(__covars[__covindex[2]])
#define seas4		(__covars[__covindex[3]])
#define seas5		(__covars[__covindex[4]])
#define seas6		(__covars[__covindex[5]])
#define S		(__x[__stateindex[0]])
#define E		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define A		(__x[__stateindex[3]])
#define R		(__x[__stateindex[4]])
#define incid		(__x[__stateindex[5]])
#define foival		(__x[__stateindex[6]])
#define Str0		(__x[__stateindex[7]])
#define Sout		(__x[__stateindex[8]])
#define Sin		(__x[__stateindex[9]])
#define cases		(__y[__obsindex[0]])
#define lik		(__lik[0])

void __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 
    if (ISNA(cases)) {
      lik = (give_log) ? 0 : 1;
    } else {
      lik = dnbinom_mu(cases, tau, rho*incid, give_log);
    }
   
}

#undef rho
#undef tau
#undef beta1
#undef beta2
#undef beta3
#undef beta4
#undef beta5
#undef beta6
#undef betat
#undef gamma
#undef sigma
#undef theta0
#undef alpha
#undef mu
#undef delta
#undef nu
#undef pop_0
#undef sig_sq
#undef S_0
#undef E_0
#undef I_0
#undef A_0
#undef R_0
#undef seas1
#undef seas2
#undef seas3
#undef seas4
#undef seas5
#undef seas6
#undef S
#undef E
#undef I
#undef A
#undef R
#undef incid
#undef foival
#undef Str0
#undef Sout
#undef Sin
#undef cases
#undef lik

/* C snippet: 'toEst' */
#define seas1		(__covars[__covindex[0]])
#define seas2		(__covars[__covindex[1]])
#define seas3		(__covars[__covindex[2]])
#define seas4		(__covars[__covindex[3]])
#define seas5		(__covars[__covindex[4]])
#define seas6		(__covars[__covindex[5]])
#define rho		(__p[__parindex[0]])
#define tau		(__p[__parindex[1]])
#define beta1		(__p[__parindex[2]])
#define beta2		(__p[__parindex[3]])
#define beta3		(__p[__parindex[4]])
#define beta4		(__p[__parindex[5]])
#define beta5		(__p[__parindex[6]])
#define beta6		(__p[__parindex[7]])
#define betat		(__p[__parindex[8]])
#define gamma		(__p[__parindex[9]])
#define sigma		(__p[__parindex[10]])
#define theta0		(__p[__parindex[11]])
#define alpha		(__p[__parindex[12]])
#define mu		(__p[__parindex[13]])
#define delta		(__p[__parindex[14]])
#define nu		(__p[__parindex[15]])
#define pop_0		(__p[__parindex[16]])
#define sig_sq		(__p[__parindex[17]])
#define S_0		(__p[__parindex[18]])
#define E_0		(__p[__parindex[19]])
#define I_0		(__p[__parindex[20]])
#define A_0		(__p[__parindex[21]])
#define R_0		(__p[__parindex[22]])
#define T_rho		(__pt[__parindex[0]])
#define T_tau		(__pt[__parindex[1]])
#define T_beta1		(__pt[__parindex[2]])
#define T_beta2		(__pt[__parindex[3]])
#define T_beta3		(__pt[__parindex[4]])
#define T_beta4		(__pt[__parindex[5]])
#define T_beta5		(__pt[__parindex[6]])
#define T_beta6		(__pt[__parindex[7]])
#define T_betat		(__pt[__parindex[8]])
#define T_gamma		(__pt[__parindex[9]])
#define T_sigma		(__pt[__parindex[10]])
#define T_theta0		(__pt[__parindex[11]])
#define T_alpha		(__pt[__parindex[12]])
#define T_mu		(__pt[__parindex[13]])
#define T_delta		(__pt[__parindex[14]])
#define T_nu		(__pt[__parindex[15]])
#define T_pop_0		(__pt[__parindex[16]])
#define T_sig_sq		(__pt[__parindex[17]])
#define T_S_0		(__pt[__parindex[18]])
#define T_E_0		(__pt[__parindex[19]])
#define T_I_0		(__pt[__parindex[20]])
#define T_A_0		(__pt[__parindex[21]])
#define T_R_0		(__pt[__parindex[22]])

void __pomp_to_trans (double *__pt, const double *__p, const int *__parindex)
{
 	T_tau = log(tau);
	T_sigma = log(sigma);
	T_gamma = log(gamma);
	T_mu = log(mu);
	T_delta = log(delta);
	T_alpha = log(alpha);
	T_sig_sq = log(sig_sq);
	T_rho = logit(rho);
	T_nu = logit(nu);
	T_theta0 = logit(theta0);
	to_log_barycentric(&T_S_0,&S_0,5); 
}

#undef seas1
#undef seas2
#undef seas3
#undef seas4
#undef seas5
#undef seas6
#undef rho
#undef tau
#undef beta1
#undef beta2
#undef beta3
#undef beta4
#undef beta5
#undef beta6
#undef betat
#undef gamma
#undef sigma
#undef theta0
#undef alpha
#undef mu
#undef delta
#undef nu
#undef pop_0
#undef sig_sq
#undef S_0
#undef E_0
#undef I_0
#undef A_0
#undef R_0
#undef T_rho
#undef T_tau
#undef T_beta1
#undef T_beta2
#undef T_beta3
#undef T_beta4
#undef T_beta5
#undef T_beta6
#undef T_betat
#undef T_gamma
#undef T_sigma
#undef T_theta0
#undef T_alpha
#undef T_mu
#undef T_delta
#undef T_nu
#undef T_pop_0
#undef T_sig_sq
#undef T_S_0
#undef T_E_0
#undef T_I_0
#undef T_A_0
#undef T_R_0

/* C snippet: 'fromEst' */
#define seas1		(__covars[__covindex[0]])
#define seas2		(__covars[__covindex[1]])
#define seas3		(__covars[__covindex[2]])
#define seas4		(__covars[__covindex[3]])
#define seas5		(__covars[__covindex[4]])
#define seas6		(__covars[__covindex[5]])
#define rho		(__p[__parindex[0]])
#define tau		(__p[__parindex[1]])
#define beta1		(__p[__parindex[2]])
#define beta2		(__p[__parindex[3]])
#define beta3		(__p[__parindex[4]])
#define beta4		(__p[__parindex[5]])
#define beta5		(__p[__parindex[6]])
#define beta6		(__p[__parindex[7]])
#define betat		(__p[__parindex[8]])
#define gamma		(__p[__parindex[9]])
#define sigma		(__p[__parindex[10]])
#define theta0		(__p[__parindex[11]])
#define alpha		(__p[__parindex[12]])
#define mu		(__p[__parindex[13]])
#define delta		(__p[__parindex[14]])
#define nu		(__p[__parindex[15]])
#define pop_0		(__p[__parindex[16]])
#define sig_sq		(__p[__parindex[17]])
#define S_0		(__p[__parindex[18]])
#define E_0		(__p[__parindex[19]])
#define I_0		(__p[__parindex[20]])
#define A_0		(__p[__parindex[21]])
#define R_0		(__p[__parindex[22]])
#define T_rho		(__pt[__parindex[0]])
#define T_tau		(__pt[__parindex[1]])
#define T_beta1		(__pt[__parindex[2]])
#define T_beta2		(__pt[__parindex[3]])
#define T_beta3		(__pt[__parindex[4]])
#define T_beta4		(__pt[__parindex[5]])
#define T_beta5		(__pt[__parindex[6]])
#define T_beta6		(__pt[__parindex[7]])
#define T_betat		(__pt[__parindex[8]])
#define T_gamma		(__pt[__parindex[9]])
#define T_sigma		(__pt[__parindex[10]])
#define T_theta0		(__pt[__parindex[11]])
#define T_alpha		(__pt[__parindex[12]])
#define T_mu		(__pt[__parindex[13]])
#define T_delta		(__pt[__parindex[14]])
#define T_nu		(__pt[__parindex[15]])
#define T_pop_0		(__pt[__parindex[16]])
#define T_sig_sq		(__pt[__parindex[17]])
#define T_S_0		(__pt[__parindex[18]])
#define T_E_0		(__pt[__parindex[19]])
#define T_I_0		(__pt[__parindex[20]])
#define T_A_0		(__pt[__parindex[21]])
#define T_R_0		(__pt[__parindex[22]])

void __pomp_from_trans (double *__p, const double *__pt, const int *__parindex)
{
 	tau = exp(T_tau);
	sigma = exp(T_sigma);
	gamma = exp(T_gamma);
	mu = exp(T_mu);
	delta = exp(T_delta);
	alpha = exp(T_alpha);
	sig_sq = exp(T_sig_sq);
	rho = expit(T_rho);
	nu = expit(T_nu);
	theta0 = expit(T_theta0);
	from_log_barycentric(&S_0,&T_S_0,5); 
}

#undef seas1
#undef seas2
#undef seas3
#undef seas4
#undef seas5
#undef seas6
#undef rho
#undef tau
#undef beta1
#undef beta2
#undef beta3
#undef beta4
#undef beta5
#undef beta6
#undef betat
#undef gamma
#undef sigma
#undef theta0
#undef alpha
#undef mu
#undef delta
#undef nu
#undef pop_0
#undef sig_sq
#undef S_0
#undef E_0
#undef I_0
#undef A_0
#undef R_0
#undef T_rho
#undef T_tau
#undef T_beta1
#undef T_beta2
#undef T_beta3
#undef T_beta4
#undef T_beta5
#undef T_beta6
#undef T_betat
#undef T_gamma
#undef T_sigma
#undef T_theta0
#undef T_alpha
#undef T_mu
#undef T_delta
#undef T_nu
#undef T_pop_0
#undef T_sig_sq
#undef T_S_0
#undef T_E_0
#undef T_I_0
#undef T_A_0
#undef T_R_0

static int __pomp_load_stack = 0;

void __pomp_load_stack_incr (void) {++__pomp_load_stack;}

void __pomp_load_stack_decr (int *val) {*val = --__pomp_load_stack;}

void R_init_pomp_6bbe27bc4528b0f995b01da8b4d49700 (DllInfo *info)
{
R_RegisterCCallable("pomp_6bbe27bc4528b0f995b01da8b4d49700", "__pomp_load_stack_incr", (DL_FUNC) __pomp_load_stack_incr);
R_RegisterCCallable("pomp_6bbe27bc4528b0f995b01da8b4d49700", "__pomp_load_stack_decr", (DL_FUNC) __pomp_load_stack_decr);
R_RegisterCCallable("pomp_6bbe27bc4528b0f995b01da8b4d49700", "__pomp_rinit", (DL_FUNC) __pomp_rinit);
R_RegisterCCallable("pomp_6bbe27bc4528b0f995b01da8b4d49700", "__pomp_stepfn", (DL_FUNC) __pomp_stepfn);
R_RegisterCCallable("pomp_6bbe27bc4528b0f995b01da8b4d49700", "__pomp_rmeasure", (DL_FUNC) __pomp_rmeasure);
R_RegisterCCallable("pomp_6bbe27bc4528b0f995b01da8b4d49700", "__pomp_dmeasure", (DL_FUNC) __pomp_dmeasure);
R_RegisterCCallable("pomp_6bbe27bc4528b0f995b01da8b4d49700", "__pomp_to_trans", (DL_FUNC) __pomp_to_trans);
R_RegisterCCallable("pomp_6bbe27bc4528b0f995b01da8b4d49700", "__pomp_from_trans", (DL_FUNC) __pomp_from_trans);
}
