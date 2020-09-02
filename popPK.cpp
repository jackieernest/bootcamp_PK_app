$PARAM
TVKA=1, TVCL = 1, TVVC=1, TVVP1=1, TVQ1=0, TVF1=1

$CMT GUT CENT P1

$MAIN
F_GUT = F1;

double KA = TVKA*exp(ETA(1));
double CL = TVCL*exp(ETA(2));
double VC = TVVC*exp(ETA(3));
double VP1 = TVVP1*exp(ETA(4));
double Q1 = TVQ1*exp(ETA(5));

double F1 = TVF1*exp(ETA(6));

double k10 = CL/VC;

double k12 = Q1/VC;
double k21 = Q1/VP1;


$OMEGA 0 0 0 0 0 0

$SIGMA @labels PROP ADD
0 0

$ODE
dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT -k12*CENT+k21*P1 - k10*CENT;
dxdt_P1 = k12*CENT-k21*P1;



$TABLE

double IPRED = CENT/VC;
double DV = IPRED*(1+PROP)+ADD;

while(DV < 0) {
  simeps();
  DV = IPRED*(1+PROP)+ADD;
}


$CAPTURE IPRED DV
