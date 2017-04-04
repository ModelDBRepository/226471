#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <iostream>
using namespace std;

//---------------CONSTANTs initialization----------------------------
//--------------- Network Geometry ------------------------------------
#define I_TC    1     //  0 - No layer, 1 - Add Layer
#define I_RE    1     //  0 - No layer, 1 - Add Layer
#define I_CX    0     //  0 - No layer, 1 - Add Layer
#define I_IN    0     //  0 - No layer, 1 - Add Layer
#define I_GB    1     //  0 - No GABAb from RE to TC, 1 - Yes
#define I_GB_RERE    0     //  0 - No GABAb from RE to RE, 1 - Yes

#define Mre       280      //  Number of LN in olfaction
#define Mtc       100      //  Number of PN in olfaction
#define Mc      1      //  Number of CX-IN cells in X direction
#define Mre1      1      //  Number of RE cells in Y direction
#define Mtc1      1      //  Number of TC cells in Y direction
#define Mc1     1      //  Number of CX-IN cells in Y direction
#define Max     ((Mre > Mtc) ? Mre : Mtc)
#define Max1    ((Mre1 > Mtc1) ? Mre1 : Mtc1)

//-------------Boundary conditions ---------------------------------------
#define BOUND    0   //  0 - flow; 1 - periodic; 9 - no boundary elements 
#define SELFING    0   //  0 - RE without self inhibition; 1 - with ... 
#define SELFINGcx  0   //  0 - CX without self excitation; 1 - with ... 

//-------------- Define the connections between cells -----------------
#define MS_RE_RE 0
#define MS_RE_RE1 0
#define MS_RE_RE_MAX  ((MS_RE_RE > MS_RE_RE1) ? MS_RE_RE : MS_RE_RE1)
#define N_RE_RE  Mre*Mre1 //(2*MS_RE_RE+1)*(2*MS_RE_RE1+1) 

#define MS_RE_TC 3
#define MS_RE_TC1 0
#define MS_RE_TC_MAX ((MS_RE_TC > MS_RE_TC1) ? MS_RE_TC : MS_RE_TC1)
#define N_RE_TC  Mre*Mre1 //(2*MS_RE_TC+1)*(2*MS_RE_TC1+1) //number of
                          //connections accepted by some TC from all RE

#define MS_TC_RE 3
#define MS_TC_RE1 0
#define MS_TC_RE_MAX ((MS_TC_RE > MS_TC_RE1) ? MS_TC_RE : MS_TC_RE1)
#define N_TC_RE  Mtc*Mtc1 //(2*MS_TC_RE+1)*(2*MS_TC_RE1+1)

#define MS_TC_TC 3
#define MS_TC_TC1 0
#define MS_TC_TC_MAX ((MS_TC_TC > MS_TC_TC1) ? MS_TC_TC : MS_TC_TC1)
#define N_TC_TC  Mtc*Mtc1 //(2*MS_TC_TC+1)*(2*MS_TC_TC1+1)

#define MS_CX_CX 0
#define MS_CX_CX1 0
#define MS_CX_CX_MAX  ((MS_CX_CX > MS_CX_CX1) ? MS_CX_CX : MS_CX_CX1)
#define N_CX_CX  (2*MS_CX_CX+1)*(2*MS_CX_CX1+1)

#define MS_CX_IN 0 
#define MS_CX_IN1 0
#define MS_CX_IN_MAX  ((MS_CX_IN > MS_CX_IN1) ? MS_CX_IN : MS_CX_IN1)
#define N_CX_IN  (2*MS_CX_IN+1)*(2*MS_CX_IN1+1)

#define MS_IN_CX 0
#define MS_IN_CX1 0
#define MS_IN_CX_MAX  ((MS_IN_CX > MS_IN_CX1) ? MS_IN_CX : MS_IN_CX1)
#define N_IN_CX  (2*MS_IN_CX+1)*(2*MS_IN_CX1+1)

#define MS_TC_CX 0
#define MS_TC_CX1 0
#define MS_TC_CX_MAX  ((MS_TC_CX > MS_TC_CX1) ? MS_TC_CX : MS_TC_CX1)
#define N_TC_CX  (2*MS_TC_CX+1)*(2*MS_TC_CX1+1)

#define MS_TC_IN 0
#define MS_TC_IN1 0
#define MS_TC_IN_MAX  ((MS_TC_IN > MS_TC_IN1) ? MS_TC_IN : MS_TC_IN1)
#define N_TC_IN  (2*MS_TC_IN+1)*(2*MS_TC_IN1+1)

#define MS_CX_TC 0 
#define MS_CX_TC1 0
#define MS_CX_TC_MAX  ((MS_CX_TC > MS_CX_TC1) ? MS_CX_TC : MS_CX_TC1)
#define N_CX_TC  (2*MS_CX_TC+1)*(2*MS_CX_TC1+1)

#define MS_CX_RE 0 
#define MS_CX_RE1 0
#define MS_CX_RE_MAX  ((MS_CX_RE > MS_CX_RE1) ? MS_CX_RE : MS_CX_RE1)
#define N_CX_RE  (2*MS_CX_RE+1)*(2*MS_CX_RE1+1)

//------------Number of ODE for earch cell -------------------------------
#define N_RE 9 //4 //7
#define N_TC 6 //12 
#define N_GB 2
#define N_GA 1

#define N_DEND   8
#define N_SOMA   3
#define N_CX     (N_DEND + N_SOMA)
#define N_IN     N_CX   //4

#define N_EQ1  (N_RE*I_RE + N_RE_RE*N_GB*I_RE*I_GB_RERE)*Mre*Mre1
#define N_EQ2  (N_TC*I_TC + N_RE_TC*N_GB*I_RE*I_TC*I_GB)*Mtc*Mtc1
#define N_EQ3  (N_CX*I_CX + N_IN*I_IN + N_IN_CX*N_GB*I_CX*I_IN)*Mc*Mc1
#define N_EQ4  (N_RE_TC*N_GA*I_RE*I_TC)*Mtc*Mtc1
#define N_EQ5  (N_RE_RE*N_GA*I_RE)*Mre*Mre1
#define N_EQ   (N_EQ1+N_EQ2+N_EQ3+N_EQ4+N_EQ5)    //  Complete number of ODE

//++++++++++++++ CURRENTs DESCRIPTION ++++++++++++++++++++++++++++++++++++++
//---------------Low-threshold Ca2+ current (RE cell)---------------------
class IT_RE {
  static double Shift, Ca_0, Cels;
  double m_inf, tau_m, h_inf, tau_h, ratio, eca, Phi_m, Phi_h, eca0;
public:
  double iT, m0, h0, Qm, Qh;
  double G_Ca;
  IT_RE(double v) {
    G_Ca = 1.75;
    Qm = 2.5; //2.5; //3; 
    Qh = 3; //2.5; //5;
    Phi_m = pow(Qm,((Cels-24)/10));
    Phi_h = pow(Qh,((Cels-24)/10));
    m0 = 1/(1 + exp(-(v + Shift + 50)/7.4));
    h0 = 1/(1 + exp((v + Shift + 78)/5));
    eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489);
    } 
  void calc(double m, double h, double &fm, double &fh, 
            double v, double cai, double x);
};

double IT_RE::Shift = 2, IT_RE::Ca_0 = 2, IT_RE::Cels = 36;

void IT_RE::calc(double m, double h, double &fm, double &fh,
                 double v, double cai, double x) {
  ratio = Ca_0/cai;
    if(ratio <= 0.)printf("\n LOG ERROR: RE: cai=%lf ratio=%lf",cai,ratio);
  eca = eca0 * log(ratio);
  iT = G_Ca*m*m*h*(v - eca);                              
  m_inf = 1/(1 + exp(-(v + 52)/7.4));
  tau_m = (3 + 1/(exp((v + 27)/10) + exp(-(v + 102)/15)))/Phi_m;
  h_inf = 1/(1 + exp((v + 80)/5));
  tau_h = (85 + 1/(exp((v + 48)/4) + exp(-(v + 407)/50)))/Phi_h;
  fm = -(m - m_inf)/tau_m;                                  
  fh = -(h - h_inf)/tau_h;                                  
}

//-------------- Ca-dependent potassium current (RE cell) --------------------
class IKCa_RE {
  static double E_KCa, Alpha, Beta, Cels;     
  double m_inf, tau_m, Tad, car;                                 
public:
  double iKCa, m0;
  double G_KCa;
  IKCa_RE(double cai) {
    G_KCa = 0.; //7; 
    Tad = pow(3,((Cels-22)/10));
    car = Alpha*cai*cai/Beta;
    m0 = car/(car + 1);  } 
  void calc(double m, double &fm, double v, double cai, double x);
};

double IKCa_RE::E_KCa = -95;
double IKCa_RE::Alpha = 48, IKCa_RE::Beta = 0.03, IKCa_RE::Cels = 36;   

void IKCa_RE::calc(double m, double &fm, double v, double cai, double x){
  iKCa = G_KCa*m*m*(v - E_KCa);                         
  car = Alpha*cai*cai/Beta;
  m_inf = car/(car + 1);
  tau_m = (1/(Beta*(car + 1)))/Tad;
  if(tau_m < 0.1) tau_m = 0.1;
  fm = -(1/tau_m)*(m - m_inf);                   
}

//-------------- Non-specific current (RE cell) --------------------
class ICAN_RE {
  static double E_CAN, Alpha, Beta, Cels;   
  double m_inf, tau_m, Tad, car;                                 
public:
  double iCAN, m0;
  double G_CAN;
  ICAN_RE(double cai) {
    G_CAN = 0; //0.06; 
    Tad = pow(3,((Cels-22)/10));
    car = Alpha*cai*cai/Beta;
    m0 = car/(car + 1);    } 
  void calc(double m, double &fm, double v, double cai, double x);
};

double ICAN_RE::E_CAN = -20;
double ICAN_RE::Alpha = 20, ICAN_RE::Beta = 0.002, ICAN_RE::Cels = 36;          

void ICAN_RE::calc(double m, double &fm, double v, double cai, double x){
  iCAN = G_CAN*m*m*(v - E_CAN);  
  car = Alpha*cai*cai/Beta;
  m_inf = car/(car + 1);
  tau_m = (1/(Beta*(car + 1)))/Tad;
  if(tau_m  < 0.1) tau_m = 0.1;
  fm = -(1/tau_m)*(m - m_inf); 
}

//--------------fast Na and K current (RE and TC cells)------------------
class INaK {
  static double Cels;
  double Alpha1, Beta1, Alpha2, Beta2, Alpha3, Beta3, v2, v2K, Phi;
  double tau_m, m_inf, tau_h, h_inf, tau_n, n_inf;
public:
  static double E_Na, E_K;
  double iK, iNa, m0, h0, n0;
  double G_Na, G_K, Vtr, VtrK;
  double S1, S2;
  INaK(double v) {
    G_K = 10;///////////////////////
    G_Na = 100;/////////////////////
    Vtr = -50;
    VtrK = -50;
    S1 = 0.32;
    S2 = 0.02;
    v2 = v - Vtr;
    v2K = v - VtrK;
    Phi = pow(3,((Cels-36)/10));
    Alpha1 = 0.32*(13 - v2)/(exp((13 - v2)/4) - 1);
    Beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5) - 1);
    m0 = Alpha1/(Alpha1 + Beta1);

    Alpha2 = 0.128*exp((17 - v2)/18);
    Beta2 = 4/(exp((40 - v2)/5) + 1);
    h0 = Alpha2/(Alpha2 + Beta2);

    Alpha3 = 0.02*(15 - v2)/(exp((15 - v2)/5) - 1);
    Beta3 = 0.5*exp((10 - v2)/40);
    n0 = Alpha3/(Alpha3 + Beta3);     } 
  void calc(double m, double h, double n, double &fm, double &fh, double &fn, 
            double v, double x);
};

double INaK::E_K = -95, INaK::E_Na = 50, INaK::Cels = 22; 

void INaK::calc(double m, double h, double n, double &fm, double &fh, double &fn,
                   double v, double x){
  v2 = v - Vtr;
  v2K = v - VtrK;
  iNa = G_Na*m*m*m*h*(v - E_Na);
  Alpha1 = S1*(13 - v2)/(exp((13 - v2)/4) - 1);
  Beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5) - 1);
  tau_m = 1/(Alpha1 + Beta1) / Phi;
  m_inf = Alpha1/(Alpha1 + Beta1);

  Alpha2 = 0.128*exp((17 - v2)/18);
  Beta2 = 4/(exp((40 - v2)/5) + 1);
  tau_h = 1/(Alpha2 + Beta2) / Phi;
  h_inf = Alpha2/(Alpha2 + Beta2);

  fm = -(m - m_inf)/tau_m;                 
  fh = -(h - h_inf)/tau_h;                 

  iK = G_K* n*n*n*n*(v - E_K);    
  Alpha3 = S2*(15 - v2K)/(exp((15 - v2K)/5) - 1);
  Beta3 = 0.5*exp((10 - v2K)/40);
  tau_n = 1/(Alpha3 + Beta3) / Phi;
  n_inf = Alpha3/(Alpha3 + Beta3);
  
  fn  = -(n - n_inf)/tau_n;                 
}

//------------------Ca-dynamics------------------------------------
class ICa {
  static double Ca_inf, K_T, K_d;
  double drive, drive0;                                 
public:
  double Taur, D;
  ICa() {Taur = 5; D = 0.1;
         drive0 = 10.0/(2.*96489.); }
  void calc(double cai, double &fcai, double iT, double x);
};

double ICa::Ca_inf = 2.4e-4;
double ICa::K_T = 0.0001, ICa::K_d = 0.0001;

void ICa::calc(double cai, double &fcai, double iT, double x) {
  drive = -drive0 * iT / D;
  if(drive < 0) drive = 0;
  fcai = drive + (Ca_inf - cai)/Taur; // - K_T*cai/(cai + K_d);
}

//------------------Low-theshold Ca2+ current (TC cell)-----------------
class IT_TC {
  static double Ca_0, Cels, Qm, Qh, Shift;
  double m_inf, tau_m, h_inf, tau_h, Phi_h, Phi_m, ratio, eca, eca0; 
//  double w, e;
public:
  double iT, m0, h0;
  double G_Ca;
  IT_TC(double v) {
     G_Ca = 2; //2;///////////////////////////////////
     Phi_m = pow(Qm,((Cels-24)/10));
     Phi_h = pow(Qh,((Cels-24)/10));
//     m0 = 1 / (1+exp(-(v+60.5)/6.2));
//     h0 = 1 / (1+exp((v+84)/4.03));
     m0 = 1 / (1+exp(-(v+59)/6.2));////////////////////////
     h0 = 1 / (1+exp((v+83)/4.0));
     eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489);  } 
  void calc(double m, double h, double &fm, double &fh,  
            double v, double cai, double x);
};

double IT_TC::Shift = 2, IT_TC::Ca_0 = 2, IT_TC::Cels = 36;
double IT_TC::Qm = 3.55, IT_TC::Qh = 3; //2.8;

void IT_TC::calc(double m, double h, double &fm, double &fh,
                 double v, double cai, double x) {
  ratio = Ca_0/cai;
    if(ratio <= 0.)printf("\n LOG ERROR: RE: cai=%lf ratio=%lf",cai,ratio);
  eca = eca0 * log(ratio);
  iT = G_Ca*m*m*h*(v - eca); 

/*    w = v * 0.001 * 2 * 96480 / (8.314*(36 + 273.16));
      if (fabs(w) > 1e-4) e = w/(exp(w)-1);
      else e = 1 - w/2;
      iT = G_Ca*m*m*h* (-0.001) * 2 * 96480 * (Ca_0 - cai*exp(w)) * e; */
/*    m_inf = 1 / (1+exp(-(v+60.5)/6.2));
      h_inf = 1 / (1+exp((v+84)/4.03));*/

  m_inf = 1 / (1+exp(-(v+59)/6.2));//////////////////////////////////////
  h_inf = 1 / (1+exp((v+83)/4.)); /////////////////////////////////////

//  iT = G_Ca*m_inf*m_inf*h*(v - eca);

  tau_m = (1/(exp(-(v+131.6)/16.7)+exp((v+16.8)/18.2)) + 0.612) / Phi_m;
//  tau_m = 0.15*m_inf*(1.7+exp(-(v+30.8)/13.5));/////////////////////
  tau_h = (30.8 + (211.4 + exp((v + Shift + 113.2)/5))/
           (1+exp((v + Shift + 84)/3.2))) / Phi_h;

  /*  if (v<-80) tau_h = exp((v+467)/66.6) / Phi_h; 
      else tau_h = (exp(-(v+21.88)/10.52)+28) / Phi_h;  */  	  

  fm = -(1/tau_m)*(m - m_inf);                                
  fh = -(1/tau_h)*(h - h_inf);
}

//----------------- h-current (TC cell) -----------------------------------
class Ih_TC {
  static double E_h, cac, pc, Shift, Cels, k2, k4, nca, nexp, taum;
  double h_inf, tau_s, alpha, beta, k3p, cc, Phi;
public:
  double ih, p10, o10, o20;
  double G_h, k1ca, ginc;

  Ih_TC(double v, double cai) {
     G_h = 0.02; //////////////
     ginc = 1.5;  //1.5;///          To decrease after-burst depolarisation
     Phi = pow(3,((Cels-36)/10));
     h_inf = 1/(1 + exp((v + 75 - Shift)/5.5));
     tau_s = (taum + 1000 / (exp((v + 71.5 - Shift)/14.2) + 
                          exp(-(v + 89 - Shift)/11.6))) / Phi;
     alpha = h_inf/tau_s;
     beta = (1 - h_inf)/tau_s;
     p10 = 1/(1 + pow((cac/cai),nca));
     o10 = 1/(1 + beta/alpha + pow((p10/pc),nexp));
     o20 = pow((p10/pc),nexp) * o10;
  }
  void calc(double o1, double p1, double o2,  
                 double &fo1, double &fp1, double &fo2,  
                 double v, double cai, double x);
};

double Ih_TC::E_h = -40, Ih_TC::cac = 0.002, Ih_TC::pc = 0.01; 
double Ih_TC::Shift = 0, Ih_TC::Cels = 36;
double Ih_TC::k2 = 0.0004 /*0.0004*/; /////To decrease after-burst depolarisation
double Ih_TC::k4 = 0.001;
double Ih_TC::nca = 4, Ih_TC::nexp = 1, Ih_TC::taum = 20;

void Ih_TC::calc(double o1, double p1, double o2, double &fo1, double &fp1, 
                 double &fo2, double v, double cai, double x) {
  ih = G_h*(o1 + ginc * o2)*(v - E_h);
  h_inf = 1/(1 + exp((v + 75 - Shift)/5.5));
  tau_s = (taum + 1000 / (exp((v + 71.5 - Shift)/14.2) + 
                          exp(-(v + 89 - Shift)/11.6))) / Phi;
  alpha = h_inf/tau_s;
  beta = (1 - h_inf)/tau_s;
  k1ca = k2 * pow((cai/cac),nca);
  k3p = k4 * pow((p1/pc),nexp);
  fo1 = alpha * (1-o1-o2) - beta * o1 + k4 * o2 - k3p * o1;
  fp1 = k1ca * (1-p1) - k2 * p1 + k4 * o2 - k3p * o1;
  fo2 = k3p * o1 - k4 * o2;
}

//----------------------Potassium A-current (TC cell)------------------------
class IA_TC {
  static double E_K, Cels;     
  double m_inf, tau_m, h_inf, tau_h, Tad;                                 
public:
  double iA, m0, h0, G_A;
  IA_TC(double v) {
    G_A = 10; //2; //////////////////////////////////////////////
    Tad = pow(3,((Cels-23.5)/10));
    m0 = 1.0 / (1+exp(-(v+60)/8.5));
    h0 = 1.0/(1+exp((v+78)/6)); } 
  void calc(double m, double h, double &fm, double &fh, double v, double x);
};

double IA_TC::Cels = 36, IA_TC::E_K = -95;

void IA_TC::calc(double m, double h, double &fm, double &fh, double v, double x){
  iA = G_A*m*m*m*m*h*(v - E_K);
  tau_m = (1.0/( exp((v+35.82)/19.69)+exp(-(v+79.69)/12.7) ) +0.37) / Tad;
  m_inf = 1.0 / (1+exp(-(v+60)/8.5));
  tau_h = 1.0/((exp((v+46.05)/5)+exp(-(v+238.4)/37.45))) / Tad;
  if(v >= -63) 
      tau_h = 19.0/Tad;
  h_inf = 1.0/(1+exp((v+78)/6));
  fm = -(1/tau_m)*(m - m_inf);                                  
  fh = -(1/tau_h)*(h - h_inf);
}

//---------------------Hight-threshold Ca2+ current (CX cell)----------------
class IHVA_CX {
  static double Shift, Ca_0, Cels, Qm, Qh, E_Ca;
  double m_inf, tau_m, h_inf, tau_h, Phi_m, Phi_h, a, b, vm;
  double ratio, eca0, eca;
public:
  double iHVA, m0, h0;
  double G_HVA;
  IHVA_CX(double v) {
    G_HVA = 0.03;
    Phi_m = pow(Qm,((Cels-23)/10));
    Phi_h = pow(Qh,((Cels-23)/10));
    vm = v + Shift;
    a = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1);
    b = 0.94*exp((-75-vm)/17);
    m0 = a/(a+b);
    a = 0.000457*exp((-13-vm)/50);
    b = 0.0065/(exp((-vm-15)/28) + 1);
    h0 = a/(a+b);
//    eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489);
    } 
  void calc(double m, double h, double &fm, double &fh, double v, double cai, double x);
};

double IHVA_CX::Shift = 0, IHVA_CX::Ca_0 = 2, IHVA_CX::E_Ca = 140; 
double IHVA_CX::Qm = 2.3, IHVA_CX::Qh = 2.3, IHVA_CX::Cels = 23; //36;

void IHVA_CX::calc(double m, double h, double &fm, double &fh,
                 double v, double cai, double x) {
//------------ECa is fixed (=140mV) instead of using Nerst eq.-------------
//  ratio = Ca_0/cai;
//  if(ratio <= 0.)printf("\n LOG ERROR: RE: cai=%lf ratio=%lf",cai,ratio);
//  eca = eca0 * log(ratio);

  iHVA = Phi_m * G_HVA * m*m*h * (v - E_Ca);
  vm = v + Shift;

  a = 0.055*(-37 - vm)/(exp((-37-vm)/3.8) - 1);
  b = 0.94*exp((-55-vm)/17);
  tau_m = (1/(a+b))/Phi_m;
  m_inf = a/(a+b);

  a = 0.000457*exp((-13-vm)/50);
  b = 0.0065/(exp((-vm-15)/28) + 1);
  tau_h = (1/(a+b))/Phi_h;
  h_inf = a/(a+b);

  fm = -(m - m_inf)/tau_m;                                  
  fh = -(h - h_inf)/tau_h;
}

//--------------Ca-dependent potassium current (CX cell)-----------------------
class IKCa_CX {
  static double E_KCa, Ra, Rb, Cels, Q, caix;     
  double m_inf, tau_m, Tad, a, b;                                 
public:
  double iKCa, m0;
  double G_KCa;
  IKCa_CX(double cai) {
    G_KCa = 0.3; 
    Tad = pow(Q,((Cels-23)/10));
    a = Ra * cai;  //------becouse caix = 1
    b = Rb;
    m0 = a/(a+b);
  }
  void calc(double m, double &fm, double v, double cai, double x);
};

double IKCa_CX::E_KCa = -90, IKCa_CX::Q = 2.3, IKCa_CX::caix = 1;
double IKCa_CX::Ra = 0.01, IKCa_CX::Rb = 0.02, IKCa_CX::Cels = 23; //36;   

void IKCa_CX::calc(double m, double &fm, double v, double cai, double x){
  iKCa = Tad * G_KCa * m * (v - E_KCa);                         

//  a = Ra * pow(cai,caix);
  a = Ra * cai;  //------becouse caix = 1
  b = Rb;
  tau_m = (1/(a+b))/Tad;
  m_inf = a/(a+b);
  fm = -(1/tau_m)*(m - m_inf);                   
}

//--------------------Potassium M-current (CX cell)-------------------------
class IKm_CX {
  static double E_Km, Ra, Rb, Cels, Q, tha, qa;     
  double m_inf, tau_m, Tad, a, b;                                 
public:
  double iKm, m0;
  double G_Km;
  IKm_CX(double v) {
    G_Km = 0.01; 
    Tad = pow(Q,((Cels-23)/10));
    a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
    b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
    m0 = a/(a+b);
  }
  void calc(double m, double &fm, double v, double x);
};

double IKm_CX::E_Km = -90, IKm_CX::Q = 2.3;
double IKm_CX::tha = -30, IKm_CX::qa = 9;
double IKm_CX::Ra = 0.001, IKm_CX::Rb = 0.001, IKm_CX::Cels = 36;   

void IKm_CX::calc(double m, double &fm, double v, double x){
  iKm = Tad * G_Km * m * (v - E_Km);
  a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
  b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
  tau_m = (1/(a+b))/Tad;
  m_inf = a/(a+b);
  fm = -(1/tau_m)*(m - m_inf);                   
}

//--------------------Fast potassium current (CX cell)-------------------
class IKv_CX {
  static double Ra, Rb, Cels, Q, tha, qa;     
  double m_inf, tau_m, Tad, a, b;                                 
public:
  static double E_Kv;
  double iKv, g_Kv, m0;
  double G_Kv;
  IKv_CX(double v) {
    G_Kv = 150; 
    Tad = pow(Q,((Cels-23)/10));

    a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
    b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
    m0 = a/(a+b);
  }
  void calc(double m, double &fm, double v, double x);
};

double IKv_CX::E_Kv = -90, IKv_CX::Q = 2.3;
double IKv_CX::tha = 25, IKv_CX::qa = 9;
double IKv_CX::Ra = 0.02, IKv_CX::Rb = 0.002, IKv_CX::Cels = 36;   

void IKv_CX::calc(double m, double &fm, double v, double x){
  g_Kv = Tad * G_Kv * m;
  iKv = g_Kv * (v - E_Kv);                         
  a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
  b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
  tau_m = (1/(a+b))/Tad;
  m_inf = a/(a+b);
  fm = -(1/tau_m)*(m - m_inf);                   
}

//----------------Fast sodium current (CX cell)--------------------------
class INa_CX {
  static double Shift, Ca_0, Cels, Qm, Qh;
  static double tha, qa, Ra, Rb, thi1, thi2, qi, thinf, qinf, Rg, Rd; 
  double m_inf, tau_m, h_inf, tau_h, Phi_m, Phi_h, a, b, vm;
  double trap0(double v, double th, double a, double q) {
        if (fabs(v/th) > 1.0e-6) {
                return ( a * (v - th) / (1 - exp(-(v - th)/q)) );
        } else {
	  return (a * q ); }
        }
public:
  static double E_Na;
  double iNa, g_Na, m0, h0;
  double G_Na;
  INa_CX(double v) {
    G_Na = 3000;
    Phi_m = pow(Qm,((Cels-23)/10));
    Phi_h = pow(Qh,((Cels-23)/10));
    vm = v + Shift;
    a = trap0(vm,tha,Ra,qa);
    b = trap0(-vm,-tha,Rb,qa);
    m0 = a/(a+b);

    a = trap0(vm,thi1,Rd,qi);
    b = trap0(-vm,-thi2,Rg,qi);
    h0 = 1/(1+exp((vm-thinf)/qinf));
    } 
  void calc(double m, double h, double &fm, double &fh, 
            double v, double x);
};

double INa_CX::Shift = -10, INa_CX::E_Na = 60; 
double INa_CX::Qm = 2.3, INa_CX::Qh = 2.3, INa_CX::Cels = 36;
double INa_CX::tha = -35, INa_CX::qa = 9;
double INa_CX::Ra = 0.182,INa_CX::Rb = 0.124; 
double INa_CX::thi1 = -50, INa_CX::thi2 = -75, INa_CX::qi = 5;
double INa_CX::thinf = -65, INa_CX::qinf = 6.2;
double INa_CX::Rg = 0.0091, INa_CX::Rd = 0.024; 

void INa_CX::calc(double m, double h, double &fm, double &fh,
                 double v, double x) {

  g_Na = Phi_m * G_Na * m*m*m*h;
  iNa = g_Na * (v - E_Na);
  vm = v + Shift;

  a = trap0(vm,tha,Ra,qa);
  b = trap0(-vm,-tha,Rb,qa);
  tau_m = (1/(a+b))/Phi_m;
  m_inf = a/(a+b);

                //"h" inactivation 
  a = trap0(vm,thi1,Rd,qi);
  b = trap0(-vm,-thi2,Rg,qi);
  tau_h = (1/(a+b))/Phi_h;
  h_inf = 1/(1+exp((vm-thinf)/qinf));

  fm = -(m - m_inf)/tau_m;                                  
  fh = -(h - h_inf)/tau_h;
}


//------------------Low-theshold Ca2+ current (LN cell)-----------------
class IT_LN {
  static double Ca_0, Cels, Qm, Qh, Shift, E_Ca;
  double m_inf, tau_m, h_inf, tau_h, Phi_h, Phi_m, ratio, eca, eca0; 
public:
  double iT, m0, h0;
  double G_Ca;
  IT_LN(double v) {
     G_Ca = 2; //2;///////////////////////////////////
     Phi_m = pow(Qm,((Cels-24)/10));
     Phi_h = pow(Qh,((Cels-24)/10));
     m0 = 1 / (1+exp(-(v+20)/6.5));////////////////////////
     h0 = 1 / (1+exp((v+25)/12.));
     eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489);  } 
  void calc(double m, double h, double &fm, double &fh,  
            double v, double cai, double x);
};

double IT_LN::Shift = 2, IT_LN::Ca_0 = 2, IT_LN::Cels = 24; //36;
double IT_LN::Qm = 3.55, IT_LN::Qh = 3, IT_LN::E_Ca = 140;

void IT_LN::calc(double m, double h, double &fm, double &fh,
                 double v, double cai, double x) {
//  ratio = Ca_0/cai;
//  if(ratio <= 0.)printf("\n LOG ERROR: RE: cai=%lf ratio=%lf",cai,ratio);
//  eca = eca0 * log(ratio);
  eca = E_Ca;
  iT = G_Ca*m*m*h*(v - eca); 

  m_inf = 1 / (1+exp(-(v+20)/6.5));//////////////////////////////////////
  h_inf = 1 / (1+exp((v+25)/12)); /////////////////////////////////////

  tau_m = 1.5; //1 + (v+30)*0.014;
    //(1/(exp(-(v+131.6)/167.)+exp((v+16.8)/182.)) + 0.612) / Phi_m;
  tau_h = 0.3*exp((v-40)/13) + 0.002*exp(-(v-60)/29);

//               (30.8 + (211.4 + exp((v + Shift + 113.2)/5))/
//               (1+exp((v + Shift + 84)/3.2))) / Phi_h;

  fm = -(1/tau_m)*(m - m_inf);                                
  fh = -(1/tau_h)*(h - h_inf);
}

//--------------fast K current (LN cells)------------------
class IKv_LN {
  static double Cels;
  double Alpha3, Beta3, v2, Phi;
  double tau_n, n_inf;
public:
  static double E_K;
  double iKv, n0;
  double G_Kv, Vtr;
  double S2;
  IKv_LN(double v) {
    G_Kv = 10;///////////////////////
    Vtr = -50;
    S2 = 0.02;
    v2 = v - Vtr;
    Phi = pow(3,((Cels-36)/10));

    Alpha3 = 0.032*(15 - v2)/(exp((15 - v2)/5) - 1);
    Beta3 = 0.5*exp((10 - v2)/40);
    n0 = Alpha3/(Alpha3 + Beta3);     } 
  void calc(double n, double &fn, double v, double x);
};

double IKv_LN::E_K = -95, IKv_LN::Cels = 22; 

void IKv_LN::calc(double n, double &fn, double v, double x){

  v2 = v - Vtr;
  iKv = G_Kv* n*n*n*n*(v - E_K);    
  Alpha3 = S2*(15 - v2)/(exp((15 - v2)/5) - 1);
  Beta3 = 0.5*exp((10 - v2)/40);
  tau_n = 1/(Alpha3 + Beta3) / Phi;
  n_inf = Alpha3/(Alpha3 + Beta3);
  
  fn  = -(n - n_inf)/tau_n;                 
}

//===================Now we'll CREATE the CELLS==================================
//-------------------RE CELL-----------------------------------------------
class RE: public IT_RE, public IKCa_RE, public ICAN_RE, 
          public INaK, public ICa {
  static double Cai0, V0;
public:
  double G_kl, G_l, E_l, S;

  RE() :IT_RE(V0), IKCa_RE(Cai0), ICAN_RE(Cai0), INaK(V0), ICa() {
        E_l = -70;
        G_l = 0.05;
        G_kl = 0.018; //0.012; //0.015;
        S = 1.43e-4;
        }
  void init(double *y) {
        y[0] = V0;
	y[1] = Cai0;
        y[2] = IT_RE::m0;
        y[3] = IT_RE::h0;
//        y[4] = IKCa_RE::m0;
//        y[5] = ICAN_RE::m0;
        y[4] = INaK::m0;
        y[5] = INaK::h0;
        y[6] = INaK::n0;
        }
  void calc(double x, double *y, double *f); 
};  

double RE::V0 = -61, RE::Cai0 = 0.0001;

void RE::calc(double x, double *y, double *f){
    IT_RE::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
//    IKCa_RE::calc(y[4], f[4], y[0],y[1], x);
//    ICAN_RE::calc(y[5], f[5], y[0],y[1], x);
//    INaK::calc(y[6], y[7], y[8], f[6], f[7], f[8], y[0], x); 
    INaK::calc(y[4], y[5], y[6], f[4], f[5], f[6], y[0], x);
    ICa::calc(y[1], f[1], iT, x);
    f[0] = -G_l * (y[0] - E_l) - iT /*- iKCa - iCAN*/ - iNa - iK
                         - G_kl * (y[0] - INaK::E_K);
}

//-------------------TC CELL-----------------------------------------------
class TC: public IT_TC, public Ih_TC, public INaK, public ICa, public IA_TC {
  static double G_l, Cai0, V0;
public:
  double G_kl, S, E_l;

  TC() :Ih_TC(V0,Cai0), IT_TC(V0), INaK(V0), ICa(), IA_TC(V0) {
        G_kl = 0.012;
        E_l = -70;
        INaK::G_Na = 90;
        INaK::G_K = 10;
        S = 2.9e-4;
        }
  void init(double *y) {
        y[0] = V0;
	y[1] = Cai0;
        y[2] = IT_TC::m0;
        y[3] = IT_TC::h0;
        y[4] = Ih_TC::o10;
        y[5] = Ih_TC::p10;
        y[6] = Ih_TC::o20;
        y[7] = INaK::m0;
        y[8] = INaK::h0;
        y[9] = INaK::n0;
        y[10] = IA_TC::m0;
        y[11] = IA_TC::h0;
        }
  void calc(double x, double *y, double *f); 
};  

double TC::Cai0 = 0.0001;
double TC::G_l = 0.01;//0.05; //0.01;///
double TC::V0 = -68;

void TC::calc(double x, double *y, double *f){
    IT_TC::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
    Ih_TC::calc(y[4], y[5], y[6], f[4], f[5], f[6], y[0], y[1], x);
    INaK::calc(y[7], y[8], y[9], f[7], f[8], f[9], y[0], x); 
    IA_TC::calc(y[10], y[11], f[10], f[11], y[0], x);
    ICa::calc(y[1], f[1], iT, x);
    f[0] = -G_l * (y[0] - E_l) - iT - ih - iNa - iK - iA
              - G_kl * (y[0] - INaK::E_K);
}

//-------------------CX CELL------------------------------------------------
//-------------------CX CELL (DENDRITE)-------------------------------------
class CX_DEND: public IHVA_CX, public IKCa_CX, public IKm_CX, 
               public INa_CX, public ICa {
  static double E_l, G_l;
public:
  double iDEND, I_Stim1;
  CX_DEND(double V0, double Cai0) :IHVA_CX(V0), IKCa_CX(Cai0), IKm_CX(V0),
                                   INa_CX(V0), ICa() { 
        INa_CX::G_Na = 1.5; 
        I_Stim1 = 0; 
        ICa::Taur = 150; //200;
        }
  void init(double V0, double Cai0, double *y) {
        y[0] = V0;
	y[1] = Cai0;
        y[2] = IHVA_CX::m0;
        y[3] = IHVA_CX::h0;
        y[4] = IKCa_CX::m0;
        y[5] = IKm_CX::m0;
        y[6] = INa_CX::m0;
        y[7] = INa_CX::h0;
  }
  void calc(double x, double *y, double *f); 
};  
double CX_DEND::E_l = -70; double CX_DEND::G_l = 1.0e3/30000;  //  mS/cm^2

void CX_DEND::calc(double x, double *y, double *f){
    IHVA_CX::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
    IKCa_CX::calc(y[4], f[4], y[0], y[1], x);
    IKm_CX::calc(y[5], f[5], y[0], x);
    INa_CX::calc(y[6], y[7], f[6], f[7], y[0], x);
    ICa::calc(y[1], f[1], iHVA, x);
    iDEND =  -G_l * (y[0] - E_l) - iHVA - iKCa - iKm - iNa + I_Stim1;     
}

//-------------------CX CELL (SOMA)-------------------------------------
class CX_SOMA: public IKv_CX, public INa_CX {
public:
  double v_SOMA, iSOMA, g1_SOMA, g2_SOMA, I_Stim2;
  CX_SOMA(double V0, double Cai0) :IKv_CX(V0), INa_CX(V0) { 
        I_Stim2 = 0; }
  void init(double V0, double Cai0, double *y) {
        v_SOMA = V0;
        y[0] = IKv_CX::m0;
        y[1] = INa_CX::m0;
        y[2] = INa_CX::h0;
        }
  void calc(double x, double *y, double *f); 
};  
void CX_SOMA::calc(double x, double *y, double *f){
    IKv_CX::calc(y[0], f[0], v_SOMA, x);
    INa_CX::calc(y[1], y[2], f[1], f[2], v_SOMA, x);
    g1_SOMA = g_Na + g_Kv;
    g2_SOMA = g_Na * INa_CX::E_Na + g_Kv * IKv_CX::E_Kv + I_Stim2; 
    iSOMA =  - iNa - iKv;     
}

//------------CX CELL (connect DENDRITE and SOMA)---------------------------
class CX: public CX_DEND, public CX_SOMA {
  static double Cai0, V0, C;
public:
  double kappa, rho, S_CX_SOMA, S_CX_DEND; 
  CX() :CX_DEND(V0,Cai0), CX_SOMA(V0,Cai0) {
        kappa = 10.0e3;      // kOm: to get mS=1/kOm 
        rho = 165; //140; //200; //165;
        }
  void init(double *y) {
        CX_DEND::init(V0, Cai0, y);
	CX_SOMA::init(V0, Cai0, y+N_DEND);
        S_CX_SOMA = 1.0e-6;  // cm^2
        S_CX_DEND = S_CX_SOMA * rho;        
        }
  void calc(double x, double *y, double *f); 
};  

double CX::Cai0 = 0.0001, CX::V0 = -68;
double CX::C = 0.75;   // uF/cm^2

void CX::calc(double x, double *y, double *f){
    CX_SOMA::calc(x, y+N_DEND, f+N_DEND);
    v_SOMA = (y[0] + kappa * S_CX_SOMA * g2_SOMA) / 
                            (1 + kappa*S_CX_SOMA * g1_SOMA);
    CX_DEND::calc(x, y, f);
    f[0] = (1.0/C) * ( iDEND + 1.0/(kappa*S_CX_DEND) * (v_SOMA - y[0]) );
}

//-------------------INTERNEURON-----------------------------------------------
class IN: public INaK, public IA_TC  {
  static double V0;
public:
  double G_kl, G_l, E_l, S, DC;

  IN() :INaK(V0), IA_TC(V0) {
        E_l = -55;
        G_l = 0.15;
        G_kl = 0.05;
        DC = 0;
        INaK::G_Na = 50;
        INaK::G_K = 10;
        INaK::Vtr = -55;
        INaK::VtrK = -45;
        S2 = 0.03;
        S = 1.43e-4;
        }
  void init(double *y) {
        y[0] = V0;
        y[1] = INaK::m0;
        y[2] = INaK::h0;
        y[3] = INaK::n0;
        y[4] = IA_TC::m0;
        y[5] = IA_TC::h0;
        }
  void calc(double x, double *y, double *f); 
};  

double IN::V0 = -61;

void IN::calc(double x, double *y, double *f){
    INaK::calc(y[1], y[2], y[3], f[1], f[2], f[3], y[0], x);
    IA_TC::calc(y[4], y[5], f[4], f[5], y[0], x);
    f[0] = -G_l * (y[0] - E_l) - iNa - iK - iA -
           G_kl * (y[0] - INaK::E_K) + DC/S;
}


//-------------------LN cell-----------------------------------------------

class LN: public INaK, public IT_LN, public ICa, public IKCa_CX, 
          public IKv_LN {
  static double V0, Cai0;
public:
  double G_kl, G_l, E_l, S, DC;

  LN() :INaK(V0), IT_LN(V0), IKCa_CX(Cai0), ICa(), IKv_LN(V0) {
    E_l = -50; //-58; //-54;
        G_l = 0.15;
        G_kl = 0.02;
        DC = 0;
        ICa::Taur = 150;
        ICa::D = 0.25;
	//        IKv_LN::Vtr = -45;
        S = 1.43e-4;
        }
  void init(double *y) {
        y[0] = V0;	
        y[1] = Cai0;
        y[2] = IT_LN::m0;
        y[3] = IT_LN::h0;
        y[4] = IKCa_CX::m0;
        y[5] = IKv_LN::n0;
        y[6] = INaK::m0;
        y[7] = INaK::h0;
        y[8] = INaK::n0;
        }
  void calc(double x, double *y, double *f); 
};  

double LN::V0 = -61, LN::Cai0 = 0.0001;

void LN::calc(double x, double *y, double *f){
    INaK::calc(y[6], y[7], y[8], f[6], f[7], f[8], y[0], x);
    IT_LN::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
    IKCa_CX::calc(y[4], f[4], y[0], y[1], x);
    ICa::calc(y[1], f[1], iT, x);
    IKv_LN::calc(y[5], f[5], y[0], x);
    f[0] = -G_l * (y[0] - E_l)  - iNa - iK - iT - iKCa - iKv -
           G_kl * (y[0] - INaK::E_K) + DC/S;
}


//---------SYNAPCES DESCRIPTION-------------------------------------------
//---------first order kiner model for GABA-A synapce---------------------
class Gaba_A {
  static double Cdur, Cmax, Deadtime, Prethresh;  
  double R, C, R0, R1;
  double lastrelease;
  double q, Rinf, Rtau;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double I, E_GABA, Alpha, Beta;
  Gaba_A() {
    E_GABA = -70; //-75; (-70 is more realistic?)
    R = 0, C = 0, R0 = 0, R1 = 0;
    Alpha = 10.5;
    Beta = 0.166;
    lastrelease = -100;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
  }
  void calc(double g_GABA_A, double x, double y_pre, double y_post);
}; 
double Gaba_A::Cdur = 0.3, Gaba_A::Cmax = 0.5, Gaba_A::Deadtime = 1;
double Gaba_A::Prethresh = -20;
void Gaba_A::calc(double g_GABA_A, double x, double y_pre, 
                    double y_post) {

        q = ((x - lastrelease) - Cdur);         
        if (q > Deadtime) {
                if (y_post > Prethresh) {        
                        C = Cmax;                
                        R0 = R;
                        lastrelease = x;
                }
        } else if (q < 0) {                     
        } else if (C == Cmax) {                  
                R1 = R;
                C = 0.;
        }
        if (C > 0) {                            
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
        }
       I = g_GABA_A * R * (y_pre - E_GABA);
}

//------second order kiner model (including G-proteins) for GABA-B synapce----
class Gaba_B {
  static double E_GABA, Cmax, Deadtime, Prethresh;   
  static double Kd, n; 
  double Gn, q;
public:
  double C, lastrelease;
  double I, r0, g0, Gn1, Cdur, K1, K2, K1K2, K3, K4, G;
  Gaba_B() {
    Cdur = 0.3; 
    K1 = 1.0; //0.52; 
//    K1K2 = 400;  //     K1K2 = K1/K2
    K2 = 0.0025; //0.0013;
    K3 = 0.1; //0.098;
    K4 = 0.06; //0.1; //0.033;
    lastrelease = -10000000;
    C = 0, r0 = 0, g0 = 0;
  }
  void calc(double r, double g, double &fr, double &fg, 
            double g_GABA_B, double x, double y_pre, double y_post);
};
double Gaba_B::E_GABA = -95, Gaba_B::Cmax = 0.5, Gaba_B::Deadtime = 1;
double Gaba_B::Prethresh = -20;
double Gaba_B::Kd = 100, Gaba_B::n = 4; 

//double Gaba_B::Prethresh = 0, Gaba_B::Kd = 100, Gaba_B::n = 4;
//double Gaba_B::K1 = 0.70; //0.63;//0.57;//0.8;//0.52;//0.09;//AAAAA 0.52; 
//double Gaba_B::K2 = 0.0035;//0.0025;//0.0022;//0.002;//0.0013;//0.0012;//AAAAA 0.0013;
//double Gaba_B::K3 = 0.15;   //0.098;//0.18;//AAAAAAAA 0.098; 
//double Gaba_B::K4 = 0.05;   //0.033;//0.034;//AAAAAAAA 0.033

void Gaba_B::calc(double r, double g, double &fr, double &fg, 
                  double g_GABA_B, double x, double y_pre, double y_post) {
        Gn = pow(g,n); 
        Gn1 = Gn/(Gn + Kd);
        G= g_GABA_B * Gn1;
        I = G * (y_pre - E_GABA);

        q = ((x - lastrelease) - Cdur);         
        if (q > Deadtime) {
                if (y_post > Prethresh) {        
                        C = Cmax;                
                        lastrelease = x; }
        } else if (q < 0) {                     
        } else if (C == Cmax) {   C = 0;  }
        fr = K1 * C * (1 - r) - r * K2;
        fg = K3 * r - K4 * g;
}

class GB: public Gaba_B {
public:
//  double y[N_GB], f[N_GB];
  GB() :Gaba_B() { }
  void init(double *y){
      lastrelease = -10000000;
      C = 0;
      y[0] = 0;
      y[1] = 0;
      }   
  void calc(double g_GABA_B, double x, double *y, double *f, 
                                             double y_pre, double y_post){
       Gaba_B::calc(y[0], y[1], f[0], f[1], g_GABA_B, x, y_pre, y_post); 
       } 
};  

//------second order kiner model (including G-proteins) for GABA-B synapce----
class Gaba_BF {
  static double E_GABA, Cmax, Deadtime, Prethresh;   
  static double Kd, n; 
  double Gn, q;
  double F, TrF, dF;
public:
  double C, lastrelease;
  double I, r0, g0, Gn1, Cdur, K1, K2, K1K2, K3, K4, G, FF;
  Gaba_BF() {
    Cdur = 0.3; 
    K1 = 1.0; //0.52; 
//    K1K2 = 400;  //     K1K2 = K1/K2
    K2 = 0.0025; //0.0013;
    K3 = 0.1; //0.098;
    K4 = 0.06; //0.1; //0.033;
    F=1;
    dF=0; //0.25;
    TrF=10000;
    lastrelease = -10000000;
    C = 0, r0 = 0, g0 = 0;
  }
  void calc(double r, double g, double &fr, double &fg, 
            double g_GABA_B, double x, double y_pre, double y_post);
};
double Gaba_BF::E_GABA = -95, Gaba_BF::Cmax = 0.5, Gaba_BF::Deadtime = 1;
double Gaba_BF::Prethresh = -20;
double Gaba_BF::Kd = 100, Gaba_BF::n = 4; 

//double Gaba_B::Prethresh = 0, Gaba_B::Kd = 100, Gaba_B::n = 4;
//double Gaba_B::K1 = 0.70; //0.63;//0.57;//0.8;//0.52;//0.09;//AAAAA 0.52; 
//double Gaba_B::K2 = 0.0035;//0.0025;//0.0022;//0.002;//0.0013;//0.0012;//AAAAA 0.0013;
//double Gaba_B::K3 = 0.15;   //0.098;//0.18;//AAAAAAAA 0.098; 
//double Gaba_B::K4 = 0.05;   //0.033;//0.034;//AAAAAAAA 0.033

void Gaba_BF::calc(double r, double g, double &fr, double &fg,  
                  double g_GABA_B, double x, double y_pre, double y_post) {
        Gn = pow(g,n); 
        Gn1 = Gn/(Gn + Kd);

        FF = g_GABA_B * F;
        G= g_GABA_B * Gn1;
        I = G * F * (y_pre - E_GABA);

        q = ((x - lastrelease) - Cdur);         
        if (q > Deadtime) {
                if (y_post > Prethresh) {        
                        C = Cmax;      
                        F = 1 + (F + dF - 1.0) * exp(-q/TrF);          
                        lastrelease = x; }
        } else if (q < 0) {                     
        } else if (C == Cmax) {   C = 0;  }
        fr = K1 * C * (1 - r) - r * K2;
        fg = K3 * r - K4 * g;
}

class GBF: public Gaba_BF {
public:
//  double y[N_GB], f[N_GB];
  GBF() :Gaba_BF() { }
  void init(double *y){
      lastrelease = -10000000;
      C = 0;
      y[0] = 0;
      y[1] = 0;
      }   
  void calc(double g_GABA_B, double x, double *y, double *f, 
                                             double y_pre, double y_post){
       Gaba_BF::calc(y[0], y[1], f[0], f[1], g_GABA_B, x, y_pre, y_post); 
       } 
};  

//------------first order kiner model for AMPA synapce---------------------
class AMPA {
  static double E_AMPA;
  static double Cdur, Cmax, Deadtime, Prethresh, Cdel; 
  static double Alpha, Beta; 
  int s;
  double R, C, R0, R1;
  double lastrelease, lastspike;
  double q, Rinf, Rtau;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double I;
  AMPA() {
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -100;
    lastspike = -100;
    s = 1;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
  }
  void calc(double g_AMPA, double x, double y_pre, double y_post);
};
double AMPA::E_AMPA = 0, AMPA::Cdur = 0.3, AMPA::Cmax = 0.5, AMPA::Deadtime = 1;
double AMPA::Cdel = 6; //synaptic delay to get some latency between depolarization
                         //of the presynaptic membrain and transmitter reliase
double AMPA::Prethresh = 0, AMPA::Alpha = 1, AMPA::Beta = 0.2;
void AMPA::calc(double g_AMPA, double x, double y_pre, 
                    double y_post) {

        q = ((x - lastrelease) - Cdur); 
        
        if(q > Deadtime) {
               if(y_post > Prethresh) {
                  if( (x - lastspike) > (Cdel + Cdur) ){
                     lastspike = x;
                     s = 1; } }  //the flag that spike was but wasn't utilized yet

               if((s == 1) && ((x - lastspike) > Cdel)) {
                  s = 0;         //spike was utilized
                  C = Cmax;                
                  R0 = R;
                  lastrelease = x;  }
        } else if (q < 0) {                     
        } else if (C == Cmax) {                  
                R1 = R;
                C = 0.;
        }

        if (C > 0) {                            
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
        }
       I = g_AMPA * R * (y_pre - E_AMPA);
}

//------------first order kiner model for AMPA synapce with Facilitation-------
class AMPA_F {
  static double E_AMPA;
  static double Cdur, Cmax, Deadtime, Prethresh, Cdel; 
  static double Alpha, Beta; 
  int s;
  double R, C, R0, R1;
  double lastrelease, lastspike;
  double q, Rinf, Rtau;
  double F, dF, TrF;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double I, FF;
  AMPA_F() {
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -100;
    lastspike = -100;
    s = 1;
    F=1;
    dF=0; //0.25;
    TrF=10000;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
  }
  void calc(double g_AMPA, double x, double y_pre, double y_post);
};
double AMPA_F::E_AMPA = 0, AMPA_F::Cdur = 0.3, AMPA_F::Cmax = 0.5, AMPA_F::Deadtime = 1;
double AMPA_F::Cdel = 6; //synaptic delay to get some latency between depolarization
                         //of the presynaptic membrain and transmitter reliase
double AMPA_F::Prethresh = 0, AMPA_F::Alpha = 1, AMPA_F::Beta = 0.2;

void AMPA_F::calc(double g_AMPA, double x, double y_pre, 
                    double y_post) {

        q = ((x - lastrelease) - Cdur); 
        
        if(q > Deadtime) {
               if(y_post > Prethresh) {
                  if( (x - lastspike) > (Cdel + Cdur) ){
                     lastspike = x;
                     s = 1; } }  //the flag that spike was but wasn't utilized yet

               if((s == 1) && ((x - lastspike) > Cdel)) {
                  s = 0;         //spike was utilized
                  C = Cmax;                
                  R0 = R;
                  F = 1; // + (F + dF - 1.0) * exptable(-q/TrF);
                  lastrelease = x;  }
        } else if (q < 0) {                     
        } else if (C == Cmax) {                  
                R1 = R;
                C = 0.;
        }

        if (C > 0) {                            
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
        }
       FF = g_AMPA * F;
       I = FF * R * (y_pre - E_AMPA);
}


//------------first order kiner model for AMPA synapce WITH depression--------------
class AMPA_D1 {
  static double E_AMPA;
  static double Cdur, Cmax, Deadtime, Prethresh, Cdel; 
  static double Alpha, Beta; 
  int s;
  double R, C, R0, R1;
  double lastrelease, lastspike, E, Use, Tr;
  double q, Rinf, Rtau;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double I;
  AMPA_D1() {
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -100;
    lastspike = -100;
    s = 1;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
    E = 1;
    Use = 0.5; 
    Tr = 500;
  }
  void calc(double g_AMPA, double x, double y_pre, double y_post);
};
double AMPA_D1::E_AMPA = 0, AMPA_D1::Cdur = 0.3, AMPA_D1::Cmax = 0.5, AMPA_D1::Deadtime = 1;
double AMPA_D1::Cdel = 0; //synaptic delay to get some latency between depolarization
                         //of the presynaptic membrain and transmitter reliase
double AMPA_D1::Prethresh = 0, AMPA_D1::Alpha = 0.94, AMPA_D1::Beta = 0.18;
void AMPA_D1::calc(double g_AMPA, double x, double y_pre, 
                    double y_post) {

        q = ((x - lastrelease) - Cdur); 
        
        if(q > Deadtime) {
               if(y_post > Prethresh) {
                  if( (x - lastspike) > (Cdel + Cdur) ){
                     lastspike = x;
                     s = 1; } }  //the flag that spike was but wasn't utilized yet

               if((s == 1) && ((x - lastspike) > Cdel)) {
                  s = 0;         //spike was utilized
                  C = Cmax;                
                  R0 = R;
                  E = 1 - (1 - E*(1-Use)) * exptable(-q/Tr);
                  lastrelease = x;  }
        } else if (q < 0) {                     
        } else if (C == Cmax) {                  
                R1 = R;
                C = 0.;
        }

        if (C > 0) {                            
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
        }
       I = g_AMPA * R * E * (y_pre - E_AMPA);
}

//---------first order kiner model for GABA-A synapce with DEPRESSION---------------------
class Gaba_A_D1 {
  static double Cdur, Cmax, Deadtime, Prethresh; 
  static double Alpha, Beta;  
  double R, C, R0, R1;
  double lastrelease, E, Use, Tr;
  double q, Rinf, Rtau;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double I, E_GABA;
  Gaba_A_D1() {
    E_GABA = -70; //-75; (-70 is more realistic?)
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -100;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
    E = 1;
    Use = 0.5; 
    Tr = 500;
  }
  void calc(double g_GABA_A, double x, double y_pre, double y_post);
}; 
double Gaba_A_D1::Cdur = 0.3, Gaba_A_D1::Cmax = 0.5, Gaba_A_D1::Deadtime = 1;
double Gaba_A_D1::Prethresh = 0, Gaba_A_D1::Alpha = 10.5, Gaba_A_D1::Beta = 0.166; 
void Gaba_A_D1::calc(double g_GABA_A, double x, double y_pre, 
                    double y_post) {

        q = ((x - lastrelease) - Cdur);         
        if (q > Deadtime) {
                if (y_post > Prethresh) {        
                        C = Cmax;                
                        R0 = R;
                        E = 1 - (1 - E*(1-Use)) * exptable(-q/Tr);
                        lastrelease = x;
                }
        } else if (q < 0) {                     
        } else if (C == Cmax) {                  
                R1 = R;
                C = 0.;
        }
        if (C > 0) {                            
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
        }
       I = g_GABA_A * E * R * (y_pre - E_GABA);
}

//-----first order kiner model for AMPA synapce used for external stimulation----
class Extern_ampa {
  static double Cdur, Cmax, Deadtime, Prethresh; 
  double R, C, R0, R1;
  double lastrelease;
  double q, Rinf, Rtau;
  double TR, w, wom, RRR;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double g, Alpha, Beta;
  Extern_ampa() {
    Alpha = 0.94;
    Beta = 0.18;
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -100;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
    TR = 10, w=0.01, wom=0;
  }
  void init(unsigned int seek) {srand(seek);}
  void calc(double g_Extern_ampa, double x);
};
double Extern_ampa::Cdur = 0.3, Extern_ampa::Cmax = 0.5, Extern_ampa::Deadtime = 1;
double Extern_ampa::Prethresh = 0;
void Extern_ampa::calc(double g_Extern_ampa, double x) {

        q = ((x - lastrelease) - Cdur);         
        if (q > Deadtime) {
                if ((x - lastrelease) > TR) {        
                        C = Cmax;                
                        R0 = R;
                        lastrelease = x;
//                        RRR = 1.0*rand()/(RAND_MAX + 1.0);
// 	  	  	  if(RRR < 0.000001) RRR = 0.000001;
//                        TR = -(log(RRR))/(w+(w/2)*sin(wom*x));
                }
        } else if (q < 0) {                     
        } else if (C == Cmax) {                  
                R1 = R;
                C = 0.;
        }
        if (C > 0) {                            
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
        }
       g = g_Extern_ampa * R;
}

//-----first order kiner model for AMPA synapce used for external stimulation----
class Extern_ampa1 {
  static double Cdur, Cmax, Deadtime, Prethresh; 
  double R, C, R0, R1;
  double lastrelease;
  double q, Rinf, Rtau;
  double TR, w, wom, RRR;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double g, Alpha, Beta, sig, Noise;
  Extern_ampa1() {
    R = 0, C = 0, R0 = 0, R1 = 0, sig = 0.8;
    Noise = 0.86;
  }
  void init(unsigned int seek) {srand(seek);}
  void calc(double g_Extern_ampa, double x);
};
void Extern_ampa1::calc(double g_Extern_ampa, double x) {
       RRR = 2.0*rand()/(RAND_MAX + 1.0) - 1.0;
       R = 1.0 + RRR/sig + 0.5*(Noise-0.86);
       g = g_Extern_ampa * R;
}

//-----first order kiner model for external pulse----
class Pulse {
  static double Cdur, Cmax, Deadtime; 
  double R, C, R0, R1, Rinf1, T1;
  double lastrelease;
  double q, Rinf, Rtau;
  double TR;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double g, Alpha, Beta, Beta1, RS;
  Pulse() {
    Alpha = 0.01;
    Beta = 0.005;
    Beta1 = 1./120.; //1./150.;
    T1 = 250; //400;
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -1500;
    TR = 2000.0;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
    Rinf1 = 1/(1+ exp(Beta1 * (-T1)));
  }
  void init() {
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
    }
  void calc(double x);
};
double Pulse::Cdur = 500, Pulse::Cmax = 1, Pulse::Deadtime = 1;
void Pulse::calc(double x) {

        q = ((x - lastrelease) - Cdur);         
        if (q > Deadtime) {
                if ( (x - lastrelease) > TR ) {        
                        C = Cmax;                
                        R0 = R;
                        lastrelease = x;
                }
        } 
        else if (q < 0) { } 
        else if (C == Cmax) {                  
                R1 = R;
                C = 0.;
        }

        if (C > 0) {                          
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
        } else {                              
         R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
//           R = (R1/Rinf1)/(1+ exp (Beta1 * (x-(lastrelease+Cdur)-T1)));
        }
     RS=R/Rinf;
}

//---------first order kiner model for gradient GABA-A synapce--------------
class Gaba_Ag {
  static double Cdur, Cmax, Deadtime, Prethresh, Sig;  
  double C;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double I, E_GABA, Alpha, Beta, G;
  Gaba_Ag() {
    E_GABA = -70; //-75; (-70 is more realistic?)
    C = 0;
    Alpha = 10;
    Beta = 0.2;
  }
  void calc(double g_GABA_A, double x, double *r, double *fr,
                               double y_pre, double y_post);
}; 
//double Gaba_Ag::Cdur = 0.3, Gaba_Ag::Cmax = 0.5, Gaba_Ag::Deadtime = 1;
double Gaba_Ag::Prethresh = -20, Gaba_Ag::Sig = 1.5;
void Gaba_Ag::calc(double g_GABA_A, double x, double *y, double *f, 
                    double y_post, double y_pre) {
   
        G = g_GABA_A * y[0];
        I = G * (y_post - E_GABA);

        C = 1/(1+exp(-(y_pre - Prethresh)/Sig));
        f[0] = Alpha * C * (1 - y[0]) - y[0] * Beta;

}

//---------first order kiner model for gradient GABA-A synapce FACILIT--------
class Gaba_AgF {
  static double Cdur, Cmax, Deadtime, Prethresh, Sig;  
  double C, F, dF, TrF, q;
  double exptable(double z)
    {
    if((z > -10) && (z < 10)) return( exp(z) );
    else return( 0 );
    }
public:
  double I, E_GABA, Alpha, Beta, FF, lastrelease, G;
  Gaba_AgF() {
    E_GABA = -70; //-75; (-70 is more realistic?)
    C = 0;
    F=1;
    dF=0; //0.2; //0.25;
    TrF=10000;
    lastrelease = -10000000;
    Alpha = 10;
    Beta = 0.2;
  }
  void calc(double g_GABA_A, double x, double *r, double *fr,
                               double y_pre, double y_post);
}; 
double Gaba_AgF::Cdur = 10, Gaba_AgF::Cmax = 0.5, Gaba_AgF::Deadtime = 1;
double Gaba_AgF::Prethresh = -20, Gaba_AgF::Sig = 1.5;
void Gaba_AgF::calc(double g_GABA_A, double x, double *y, double *f, 
                    double y_post, double y_pre) {

        FF = g_GABA_A * F;
        G = FF * y[0];
        I = G * (y_post - E_GABA);

        C = 1/(1+exp(-(y_pre - Prethresh)/Sig));
        q = ((x - lastrelease) - Cdur);         
        if (q > Deadtime) {
                if (y_pre > Prethresh) {        
                        F = 1 + (F + dF - 1.0) * exptable(-q/TrF);
                        lastrelease = x;
                        }
                } 
        f[0] = Alpha * C * (1 - y[0]) - y[0] * Beta;
}


//------second order kiner model (including G-proteins) for GABA-B synapce---
class Gaba_Bg {
  static double E_GABA, Cmax, Deadtime, Prethresh, Sig;   
  static double Kd, n; 
  double Gn, q;
public:
  double C, lastrelease;
  double I, r0, g0, Gn1, Cdur, K1, K2, K1K2, K3, K4;
  Gaba_Bg() {
    Cdur = 0.3; 
    K1 = 0.52; 
    K2 = 0.0013; 
    K3 = 0.098;
    K4 = 0.033;
    lastrelease = -10000000;
    C = 0, r0 = 0, g0 = 0;
  }
  void calc(double r, double g, double &fr, double &fg, 
            double g_GABA_B, double x, double y_pre, double y_post);
};
double Gaba_Bg::E_GABA = -95, Gaba_Bg::Cmax = 0.5, Gaba_Bg::Deadtime = 1;
double Gaba_Bg::Prethresh = -20,  Gaba_Bg::Sig = 0.1;
double Gaba_Bg::Kd = 100, Gaba_Bg::n = 4; 

void Gaba_Bg::calc(double r, double g, double &fr, double &fg, 
                  double g_GABA_B, double x, double y_pre, double y_post) {
        Gn = pow(g,n); 
        Gn1 = Gn/(Gn + Kd);
        I = g_GABA_B * Gn1 * (y_pre - E_GABA);

        C = 1/(1+exp(-(y_post - Prethresh)/Sig));
       
        fr = K1 * C * (1 - r) - r * K2;
        fg = K3 * r - K4 * g;
}

class GBg: public Gaba_Bg {
public:
//  double y[N_GB], f[N_GB];
  GBg() :Gaba_Bg() { }
  void init(double *y){
      lastrelease = -10000000;
      C = 0;
      y[0] = 0;
      y[1] = 0;
      }   
  void calc(double g_GABA_B, double x, double *y, double *f, 
                                             double y_pre, double y_post){
       Gaba_Bg::calc(y[0], y[1], f[0], f[1], g_GABA_B, x, y_pre, y_post); 
       } 
};  

//+++++++++++++++++++ MAIN PROGRAM +++++++++++++++++++++++++++++++++++++++++++

//----------external functuons----------------------------------------------
void rk(unsigned, void (unsigned, double, double*, double*), double, double, 
                                double*, double*, double*, double*);
void fun(unsigned, double, double*, double*);
void solout (long, double, double, double*, unsigned, int*);
//----------external variables---------------------------------------------
double g_gaba_a, g_gaba_a1, g_gaba_b, g_gaba_b_re_re, g_ampa, g_ampa1;
double g_ampa_cx_cx,  g_ampa_cx_tc, g_ampa_cx_re, g_ampa_tc_cx; 
double g_ampa_cx_in,  g_gaba_a_in_cx, g_gaba_b_in_cx, g_ampa_tc_in; 
double g_ext_tc, g_ext_re, g_ext_cx, g_ext_in;
double g_ext1, g_ext2, g_ext3, g_ext4, Noise, tt;
double SynapsesTCRE[Mre*Mre1*N_TC_RE][5], SynapsesTCTC[Mtc*Mtc1*N_TC_TC][5], SynapsesRETC[Mtc*Mtc1*N_RE_TC][5],  SynapsesRERE[Mre*Mre1*N_RE_RE][5], SynapsesRETCB[Mtc*Mtc1*N_RE_TC][5]; 

FILE *f1, *f2, *f3, *f4, *f6, *f7, *f8, *f9, *f10, *f11, *f12, *f13, *f22, *f222, *fNoise, *fSYNre, *fSYNtc, *fsRERE, *fsRETC, *fsTCRE, *fsTCTC, *fsRETCB, *fre, *ftc,*fC_PN,*fC_LN,*fC_LNLN,*fC_LNPN,*fC_PNLN,*fC_PNPN, *ftest1, *ftest2, *fStimPN, *fStimLN,*fGABA_A_TCIP, *fGABA_B_TCIP,*fAMPA_TCIP;

double xx,xx1,xx2,xx3,xx4,xx5;

int no_re[Mre][Mre1][N_RE], no_tc[Mtc][Mtc1][N_TC];
int no_g[Mtc][Mtc1][N_RE_TC][N_GB], no_gb_rere[Mre][Mre1][N_RE_RE][N_GB];
int no_cx[Mc][Mc1][N_CX], no_in[Mc][Mc1][N_IN], no_gcx[Mc][Mc1][N_IN_CX][N_GB];
int no_ga_retc[Mtc][Mtc1][N_RE_TC][N_GA], no_ga_rere[Mre][Mre1][N_RE_RE][N_GA];
double C_RERE[Mre][Mre1][Mre][Mre1], C_RETC[Mre][Mre1][Mtc][Mtc1];
double C_TCRE[Mtc][Mtc1][Mre][Mre1], C_TCTC[Mtc][Mtc1][Mtc][Mtc1];
double C_RE[Mre][Mre1], C_TC[Mtc][Mtc1], C_RE1[Mre][Mre1], C_TC1[Mtc][Mtc1];
int k_re_re[Mre][Mre1], k_re_tc[Mtc][Mtc1], k_tc_re[Mre][Mre1], k_tc_tc[Mtc][Mtc1];
double inh_a[Mtc][Mtc1],inh_b[Mtc][Mtc1],exc_AMPA[Mtc][Mtc1];
//int **C_RE, **C_TC;
int KMAXrere, KMAXretc, KMAXtcre, KMAXtctc, KMAXretcB;
/*
//----------external variables ---------------------------------------------
double g_gaba_a, g_gaba_a1, g_gaba_b, g_gaba_b_re_re, g_ampa, g_ampa1;
double g_ampa_cx_cx,  g_ampa_cx_tc, g_ampa_cx_re, g_ampa_tc_cx; 
double g_ampa_cx_in,  g_gaba_a_in_cx, g_gaba_b_in_cx, g_ampa_tc_in; 
double g_ext_tc, g_ext_re, g_ext_cx, g_ext_in;
double g_ext1, g_ext2, g_ext3, g_ext4, Noise, tt;
double SynapsesTCRE[Mre*Mre1*N_TC_RE][5], SynapsesTCTC[Mtc*Mtc1*N_TC_TC][5], SynapsesRETC[Mtc*Mtc1*N_RE_TC][5],  SynapsesRERE[Mre*Mre1*N_RE_RE][5], SynapsesRETCB[Mtc*Mtc1*N_RE_TC][5]; 

FILE *f1, *f2, *f3, *f4, *f6, *f7, *f8, *f9, *f10, *f11, *f12, *f13, *f22, *f222, *fNoise, *fSYNre, *fSYNtc, *fsRERE, *fsRETC, *fsTCRE, *fsTCTC, *fsRETCB, *fre, *ftc;

int no_re[Mre][Mre1][N_RE], no_tc[Mtc][Mtc1][N_TC];
int no_g[Mtc][Mtc1][N_RE_TC][N_GB], no_gb_rere[Mre][Mre1][N_RE_RE][N_GB];
int no_cx[Mc][Mc1][N_CX], no_in[Mc][Mc1][N_IN], no_gcx[Mc][Mc1][N_IN_CX][N_GB];
int no_ga_retc[Mtc][Mtc1][N_RE_TC][N_GA], no_ga_rere[Mre][Mre1][N_RE_RE][N_GA];
long int C_RERE[Mre][Mre1][Mre][Mre1], C_RETC[Mre][Mre1][Mtc][Mtc1];
long int C_TCRE[Mtc][Mtc1][Mre][Mre1], C_TCTC[Mtc][Mtc1][Mtc][Mtc1];
long int C_RE[Mre][Mre1], C_TC[Mtc][Mtc1], C_RE1[Mre][Mre1], C_TC1[Mtc][Mtc1];
long int k_re_re[Mre][Mre1], k_re_tc[Mtc][Mtc1], k_tc_re[Mre][Mre1], k_tc_tc[Mtc][Mtc1];
//int **C_RE, **C_TC;
long int KMAXrere, KMAXretc, KMAXtcre, KMAXtctc, KMAXretcB;
*/
//----------external classes (beginning of initialization)------------------
Gaba_AgF      *g_re_re[Mre][Mre1];
Gaba_AgF      *g_re_tc[Mtc][Mtc1];
GBF          *gb[Mtc][Mtc1];
GBF          *gb_re_re[Mre][Mre1];
AMPA_F         *a_tc_re[Mre][Mre1];
AMPA_F         *a_tc_tc[Mtc][Mtc1];
LN           re_cell[Mre][Mre1];
IN           tc_cell[Mtc][Mtc1];

AMPA         *a_cx_cx[Mc][Mc1];
AMPA         *a_cx_in[Mc][Mc1];
Gaba_A       *ga_in_cx[Mc][Mc1];
GB           *gb_in_cx[Mc][Mc1];
CX           cx_cell[Mc][Mc1];
CX           in_cell[Mc][Mc1];         

AMPA         *a_cx_tc[Mtc][Mtc1];
AMPA         *a_cx_re[Mre][Mre1];
AMPA         *a_tc_cx[Mc][Mc1];
AMPA         *a_tc_in[Mc][Mc1];

Extern_ampa1  a_ext1, a_ext2, a_ext3, a_ext4;
Pulse aP, aPRE[Mre][Mre1];

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc,char **argv)
{
//---------allocate place for ALL variables and functions-----------
  double y_ini[N_EQ], f_ini[N_EQ];

//---------allocate place for TWO temporal arrays using by RK.C solver 
  double y1[N_EQ], y2[N_EQ];

//---------general parameters----------------------------------------------
  printf("\n ... reading input parameters");
  double t = 0, tmax, tmax1, t3D, ttime, TAU, *Noi;
  double h = 0.04, red, t_ext, g_Extern_ampa3, g_Extern_ampa4, R, r[2];
  int i, j, k, l, i1, j1, ii = 0, i_inp, i_gip = 0, ih, ni;
  double g_GABA_A, g_GABA_A1, g_GABA_B, g_AMPA, g_AMPA1, g_GABA_B_RE_RE;
  double g_AMPA_CX_CX, g_AMPA_CX_TC, g_AMPA_CX_RE, g_AMPA_TC_CX;
  double g_AMPA_CX_IN, g_GABA_A_IN_CX, g_GABA_B_IN_CX, g_AMPA_TC_IN;
  double g_Extern_ampa1, g_Extern_ampa2, g_Extern_ampa, av=0, avr=0;
  int k1=0, p=0, pmax=0, nt=17501;
  int *vm[Mtc][Mtc1];
  double tsp[Mtc][Mtc1], *v[Mtc][Mtc1], *vs[Mtc][Mtc1];
  unsigned ic=10;
  int tc;
//---------solver parameters (DOPRI)-------------------------------------------
  int      res, iout = 2, itoler = 0; 
  double   atoler = 1.0E-5, rtoler = 1.0E-5; 

//----------arrays initialization----------------------------------------------
  for(i=0; i<N_EQ; i++){
    y_ini[i] = 0, f_ini[i] = 0; 
    y1[i] = 0, y2[i] = 0; }

//-------parameter initialization (from file)----------------------------------
if (argc <= 1) {
  puts("Command parameters");
  puts("-----------------------");
  puts("Input File"); }

if (!(f1=fopen(argv[1],"r"))) {
   printf("%s doesn't exist\n",argv[1]);
   exit(0); }

  fscanf(f1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", 
      &tmax, &t3D, &ttime, &t_ext, 
      &g_GABA_A, &g_GABA_A1, &g_GABA_B, &g_GABA_B_RE_RE, &g_AMPA, &g_AMPA1,
      &g_AMPA_CX_CX, &g_AMPA_CX_TC, &g_AMPA_CX_RE, &g_AMPA_TC_CX, 
      &g_AMPA_CX_IN, &g_GABA_A_IN_CX, &g_GABA_B_IN_CX, &g_AMPA_TC_IN,  
      &g_Extern_ampa1, &g_Extern_ampa2, &g_Extern_ampa3, &g_Extern_ampa4, 
      &i_inp);
  printf("\n param: %lf %lf %lf %lf A=%lf A1=%lf B=%lf B1=%lf AM=%lf AM1=%lf AM_CX_CX=%lf AM_CX_TC=%lf AM_CX_RE=%lf AM_TC_CX=%lf AM_CX_IN=%lf GA_IN_CX=%lf GB_IN_CX=%lf AM_TC_IN=%lf  %lf %lf %lf %lf %d", 
       tmax, t3D, ttime, t_ext, 
       g_GABA_A, g_GABA_A1, g_GABA_B, g_GABA_B_RE_RE, g_AMPA, g_AMPA1, 
       g_AMPA_CX_CX, g_AMPA_CX_TC, g_AMPA_CX_RE, g_AMPA_TC_CX,
       g_AMPA_CX_IN, g_GABA_A_IN_CX, g_GABA_B_IN_CX, g_AMPA_TC_IN,  
       g_Extern_ampa1, g_Extern_ampa2, g_Extern_ampa3, g_Extern_ampa4, 
       i_inp);
  fclose(f1);

  cout << "\n RE-RE_MAX=" << MS_RE_RE_MAX << " RE-TC_MAX=" << MS_RE_TC_MAX << 
                                          " TC-RE_MAX=" << MS_TC_RE_MAX;

//----------classes initialization (continue)----------------------------
  for(i=0; i < Mre; ++i)
  for(j=0; j < Mre1; ++j){
    g_re_re[i][j] = new Gaba_AgF[N_RE_RE];
    a_tc_re[i][j] = new AMPA_F[N_TC_RE];
    a_cx_re[i][j] = new AMPA[N_CX_RE]; 
    if(I_GB_RERE == 1) gb_re_re[i][j] = new GBF[N_RE_RE];
  }

  for(i=0; i < Mtc; ++i)
  for(j=0; j < Mtc1; ++j){
    g_re_tc[i][j] = new Gaba_AgF[N_RE_TC];
    if(I_GB == 1) gb[i][j] = new GBF[N_RE_TC];
    a_tc_tc[i][j] = new AMPA_F[N_TC_TC];
    a_cx_tc[i][j] = new AMPA[N_CX_TC];
  }

  for(i=0; i < Mc; ++i)
  for(j=0; j < Mc1; ++j){ 
    if(I_CX == 1) a_cx_cx[i][j] = new AMPA[N_CX_CX];
    if(I_CX == 1 && I_IN == 1) a_cx_in[i][j] = new AMPA[N_CX_IN];
    if(I_CX == 1 && I_IN == 1) ga_in_cx[i][j] = new Gaba_A[N_IN_CX];
    if(I_CX == 1 && I_IN == 1) gb_in_cx[i][j] = new GB[N_IN_CX];

    if(I_CX == 1 && I_TC == 1) a_tc_cx[i][j] = new AMPA[N_TC_CX];
    if(I_IN == 1 && I_TC == 1) a_tc_in[i][j] = new AMPA[N_TC_IN];
  }

   for(i=0; i < Mc; ++i)
   for(j=0; j < Mc1; ++j)
     in_cell[i][j].rho = 50;         //

   for(i=0; i < Mtc; ++i) 
      for(j=0; j < Mtc1; ++j){
            v[i][j] = new double[nt];
            vs[i][j] = new double[nt];
            vm[i][j] = new int[nt];
}
   for(i=0; i < Mtc; ++i) 
      for(j=0; j < Mtc1; ++j)
         for(k=0; k < nt; ++k){
            v[i][j][k] = 0.0;
            vm[i][j][k] = 0;
            vs[i][j][k] = 0.0;
         }

   Noi = new double[2500000];

/* for(i=0; i < Mre; ++i) 
      for(j=0; j < Mre1; ++j){
         for(k=0; k < Mre; ++k) 
            C_RERE[i][j][k] = new int[Mre1];
         for(k=0; k < Mtc; ++k) 
            C_RETC[i][j][k] = new int[Mtc1];
}
   for(i=0; i < Mtc; ++i) 
      for(j=0; j < Mtc1; ++j){
         for(k=0; k < Mre; ++k) 
            C_TCRE[i][j][k] = new int[Mre1];
         for(k=0; k < Mtc; ++k) 
            C_TCTC[i][j][k] = new int[Mtc1];
}
   C_TC = (int**) new int[Mtc];
   for(i=0; i < Mtc; ++i) 
      C_TC[i] = new int[Mtc1];
   for(i=0; i < Mtc; ++i) 
      for(j=0; j < Mtc1; ++j)
          printf("\n TC %d %d %d", i,j,C_TC[i][j]);

   C_RE = (int**) new int[Mre];
   for(i=0; i < Mre; ++i) 
      C_RE[i] = new int[Mre1];
   for(i=0; i < Mre; ++i) 
      for(j=0; j < Mre1; ++j)
          printf("\n RE %d %d %d", i,j,C_RE[i][j]);
*/

//----------creating the integer arrays containing the addresses--------
//----------of  ALL internal variables for ALL objects RE, TC ----------
//----------and GB classes (e.g., no_re[i][j][k] is the address --------
//----------of the variable y[k] for the object re_cell[i][j]) ---------
//----------NOTE: this is the relative adresses and you should ---------
//----------add the REAL address of the first element of the -----------
//----------original 1D array to use these arrays-----------------------
  for(i=0; i < Mre; ++i)
  for(j=0; j < Mre1; ++j)
  for(k=0; k < N_RE; ++k)
     no_re[i][j][k] = k + (j + i*Mre1) * N_RE;
  for(i=0; i < Mtc; ++i)
  for(j=0; j < Mtc1; ++j)
  for(k=0; k < N_TC; ++k)
     no_tc[i][j][k] = Mre*Mre1*N_RE*I_RE + k + (j + i*Mtc1) * N_TC;  
  for(i=0; i < Mtc; ++i)
  for(j=0; j < Mtc1; ++j)
  for(k=0; k < N_RE_TC; ++k)
  for(l=0; l < N_GB; ++l)
     no_g[i][j][k][l] = Mre*Mre1*N_RE*I_RE + Mtc*Mtc1*N_TC*I_TC + l + 
                                   (k + (j + i*Mtc1)*N_RE_TC) * N_GB;
  for(i=0; i < Mtc; ++i)
  for(j=0; j < Mtc1; ++j)
  for(k=0; k < N_RE_TC; ++k)
  for(l=0; l < N_GA; ++l)
     no_ga_retc[i][j][k][l] = Mre*Mre1*N_RE*I_RE + Mtc*Mtc1*N_TC*I_TC + 
                              Mtc*Mtc1*N_RE_TC*N_GB*I_RE*I_TC*I_GB +    
                              l + (k + (j + i*Mtc1)*N_RE_TC) * N_GA;
  for(i=0; i < Mre; ++i)
  for(j=0; j < Mre1; ++j)
  for(k=0; k < N_RE_RE; ++k)
  for(l=0; l < N_GA; ++l)
     no_ga_rere[i][j][k][l] = Mre*Mre1*N_RE*I_RE + Mtc*Mtc1*N_TC*I_TC + 
                              Mtc*Mtc1*N_RE_TC*N_GB*I_RE*I_TC*I_GB + 
                              Mtc*Mtc1*N_RE_TC*N_GA*I_RE*I_TC +
                              l + (k + (j + i*Mre1)*N_RE_RE) * N_GA;
  for(i=0; i < Mre; ++i)
  for(j=0; j < Mre1; ++j)
  for(k=0; k < N_RE_RE; ++k)
  for(l=0; l < N_GB; ++l)
     no_gb_rere[i][j][k][l] = Mre*Mre1*N_RE*I_RE + Mtc*Mtc1*N_TC*I_TC +
                              Mtc*Mtc1*N_RE_TC*N_GB*I_RE*I_TC*I_GB + 
                              Mtc*Mtc1*N_RE_TC*N_GA*I_RE*I_TC +
                              Mre*Mre1*N_RE_RE*N_GA*I_RE +
                              + l + (k + (j + i*Mre1)*N_RE_RE) * N_GB;

  for(i=0; i < Mc; ++i)
  for(j=0; j < Mc1; ++j)
  for(k=0; k < N_CX; ++k)
     no_cx[i][j][k] = Mre*Mre1*N_RE*I_RE + Mtc*Mtc1*N_TC*I_TC + 
                              Mtc*Mtc1*N_RE_TC*N_GB*I_RE*I_TC*I_GB +
                              Mtc*Mtc1*N_RE_TC*N_GA*I_RE*I_TC + 
                              Mre*Mre1*N_RE_RE*N_GA*I_RE +
                              Mre*Mre1*N_RE_RE*N_GB*I_RE*I_GB_RERE +
                                    k + (j + i*Mc1) * N_CX;
  for(i=0; i < Mc; ++i)
  for(j=0; j < Mc1; ++j)
  for(k=0; k < N_IN; ++k)
     no_in[i][j][k] = Mre*Mre1*N_RE*I_RE + Mtc*Mtc1*N_TC*I_TC + 
                               Mtc*Mtc1*N_RE_TC*N_GB*I_RE*I_TC*I_GB +
                               Mtc*Mtc1*N_RE_TC*N_GA*I_RE*I_TC + 
                               Mre*Mre1*N_RE_RE*N_GA*I_RE +
                               Mre*Mre1*N_RE_RE*N_GB*I_RE*I_GB_RERE +
                               Mc*Mc1*N_CX*I_CX + k + (j + i*Mc1) * N_IN;
  for(i=0; i < Mc; ++i)
  for(j=0; j < Mc1; ++j)
  for(k=0; k < N_IN_CX; ++k)
  for(l=0; l < N_GB; ++l)
     no_gcx[i][j][k][l] = Mre*Mre1*N_RE*I_RE + Mtc*Mtc1*N_TC*I_TC + 
                          Mtc*Mtc1*N_RE_TC*N_GB*I_RE*I_TC*I_GB +
                          Mtc*Mtc1*N_RE_TC*N_GA*I_RE*I_TC + 
                          Mre*Mre1*N_RE_RE*N_GA*I_RE +
                          Mre*Mre1*N_RE_RE*N_GB*I_RE*I_GB_RERE +
                          Mc*Mc1*N_CX*I_CX + Mc*Mc1*N_IN*I_IN +
                          l + (k + (j + i*Mc1)*N_IN_CX) * N_GB;
                                   
//---variable initialization (additional for standart constructor)----------
  for(i=0; i<Mre; i++)
    for(j=0; j<Mre1; j++){
      if(I_RE == 1) re_cell[i][j].init(y_ini+no_re[i][j][0]);
    }

  for(i=0; i<Mtc; i++)
    for(j=0; j<Mtc1; j++){
      if(I_TC == 1) tc_cell[i][j].init(y_ini+no_tc[i][j][0]);
//      if(I_RE == 1 && I_TC == 1 && I_GB == 1) for(k=0; k < N_RE_TC; k++){
//                           gb[i][j][k].init(y_ini+no_g[i][j][k][0]); }    
    }

  for(i=0; i<Mc; i++)
    for(j=0; j<Mc1; j++){
      if(I_CX == 1) cx_cell[i][j].init(y_ini+no_cx[i][j][0]); 
      if(I_IN == 1) in_cell[i][j].init(y_ini+no_in[i][j][0]);
      if(I_CX == 1 && I_IN == 1) for(k=0; k < N_IN_CX; k++){
                       gb_in_cx[i][j][k].init(y_ini+no_gcx[i][j][k][0]); }     
    }
  a_ext1.init(8);
  a_ext2.init(9);
  a_ext3.init(10);
  a_ext4.init(11);
    for(j = 0; j < Mtc1; ++j)
       for(i = 0; i < Mtc; ++i) 
          tsp[i][j] = -100; 


printf("\n NOISE");
fNoise = fopen("time_Hz100Syn200", "r");
i=0;
    while(tt< 99000.){
       fscanf(fNoise,"%lf %lf",&tt, &Noise);
       Noi[i] = Noise;
       //       printf("\n %lf %lf",tt,Noise);
       i=i+1;
    }
fclose(fNoise);
//----------here we are changing some variables----------------------


  for(i=0; i < Mre; ++i)
  for(j=0; j < Mre1; ++j){
     re_cell[i][j].G_kl = 0.01; //0.0015; //0.005;
  }
  for(i=0; i < Mtc; ++i)
  for(j=0; j < Mtc1; ++j){
    tc_cell[i][j].G_kl = 0.04; //0.01; //0.012; //0.015;
  }
  for(i=0; i < Mre; ++i)
  for(j=0; j < Mre1; ++j){
     re_cell[i][j].G_KCa = 0.25;
     //     re_cell[i][j].G_Ca = 0;
     re_cell[i][j].G_Kv = 7;
  }
  a_ext1.sig = 1000000000.0;
  a_ext2.sig = 1000000000.0;

  //  a_ext3.sig = 1000000000.0;
  //  a_ext4.sig = 1000000000.0;
//---changes of the parameters to get variability------------------------
     r[0] = 2.0 * rand()/(RAND_MAX + 1.0) - 1.0;
     r[1] = 2.0 * rand()/(RAND_MAX + 1.0) - 1.0;

srand(1);
  if(I_RE == 1)
  for(i=0; i < Mre; ++i)
   for(j=0; j < Mre1; ++j){
A:   R = 2.0 * rand()/(RAND_MAX + 1.0) - 1.0;
     if((R < -1) || (R > 1)) cout << "\n CHECK RANDOM !!!!!!!!";
//     if(R < -0.6) goto A;
//     if((r[0]*r[1]>0) && (r[1]*R>0)) goto A;
     while((r[0]*r[1]>0) && (r[1]*R>0)) 
                         {R = 2.0 * rand()/(RAND_MAX + 1.0) - 1.0;}
     r[1]=r[0]; r[0]=R;
     re_cell[i][j].G_kl = re_cell[i][j].G_kl + R * 0.005; //0.015; //0.02; 
     cout << "\n R=" << R << "  G_l(RE)=" << re_cell[i][j].G_kl; }
cout << "\n -------------------------------------------------";

srand(3);
  if(I_TC == 1) 
  for(i=0; i < Mtc; ++i)
   for(j=0; j < Mtc1; ++j){
     R = 2.0 * rand()/(RAND_MAX + 1.0) - 1.0;
     if((R < -1) || (R > 1)) cout << "\n CHECK RANDOM !!!!!!!!";
     while((r[0]*r[1]>0) && (r[1]*R>0)) {R = 2.0 * rand()/(RAND_MAX + 1.0) - 1.0;}
     r[1]=r[0]; r[0]=R;
     tc_cell[i][j].G_kl = tc_cell[i][j].G_kl + R * 0.005; 
     cout << "\n R=" << R << "  G_kl(TC)=" << tc_cell[i][j].G_kl; }
cout << "\n -------------------------------------------------";

/* srand(2);
   if(I_TC == 1)
   for(i=0; i < Mtc; ++i)
    for(j=0; j < Mtc1; ++j){
     R = 2.0 * rand()/(RAND_MAX + 1.0) - 1.0;
     if((R < -1) || (R > 1)) cout << "\n CHECK RANDOM !!!!!!!!";
     while((r[0]*r[1]>0) && (r[1]*R>0)) {R = 2.0 * rand()/(RAND_MAX + 1.0) - 1.0;}
     r[1]=r[0]; r[0]=R;
     tc_cell[i][j].G_h = tc_cell[i][j].G_h + R * 0.002; 
     cout << "\n R=" << R << "  G_h(TC)=" << tc_cell[i][j].G_h; }
cout << "\n -------------------------------------------------";
*/

for(i=0; i<Mtc; i++)
for(j=0; j<Mtc1; j++){
  R = 1.0 * rand()/(RAND_MAX + 1.0);
  if(R < 0.1) { tc_cell[i][j].DC=0.0001;
                printf("\n TC %d %d", i,j);}
  }

for(i=0; i<Mre; i++)
for(j=0; j<Mre1; j++){
  R = 1.0 * rand()/(RAND_MAX + 1.0);
  if(R < 0.1) { re_cell[i][j].DC=0.0001;
                printf("\n RE %d %d", i,j);}
  }

//--------------end variability------------------------------------
printf("\n Zero connections");
srand(7);
for(i=0; i<Mtc; i++)
for(j=0; j<Mtc1; j++) {
   printf("\n %d %d %d", i,j,C_TC[i][j]);
   C_TC[i][j]=0;
}

for(i=0; i<Mre; i++)
for(j=0; j<Mre1; j++) {
   printf("\n %d %d %d", i,j,C_RE[i][j]);
   C_RE[i][j]=0;
}

for(i=0; i<Mre; i++)
for(j=0; j<Mre1; j++)
for(i1=0; i1<Mre; i1++)
for(j1=0; j1<Mre1; j1++) C_RERE[i][j][i1][j1]=0;

for(i=0; i<Mre; i++)
for(j=0; j<Mre1; j++)
for(i1=0; i1<Mtc; i1++)
for(j1=0; j1<Mtc1; j1++) C_RETC[i][j][i1][j1]=0;

for(i=0; i<Mtc; i++)
for(j=0; j<Mtc1; j++)
for(i1=0; i1<Mre; i1++)
for(j1=0; j1<Mre1; j1++) C_TCRE[i][j][i1][j1]=0;

for(i=0; i<Mtc; i++)
for(j=0; j<Mtc1; j++)
for(i1=0; i1<Mtc; i1++)
for(j1=0; j1<Mtc1; j1++) C_TCTC[i][j][i1][j1]=0;

printf("\n Start connections \n");

fC_PN = fopen("C_PN.txt","r");

for(i=0; i<Mtc; i++)
for(j=0; j<Mtc1; j++){
   fscanf(fC_PN,"%lf",&xx);
   C_TC[i][j] = xx;
// printf("\n C_TC=%lf",C_TC[i][j]);
    printf("%lf\n",C_TC[i][j]);
  }
 fclose(fC_PN);

fC_LN = fopen("C_LN.txt","r");

for(i=0; i<Mre; i++)
for(j=0; j<Mre1; j++){
   fscanf(fC_LN,"%lf",&xx5);
   C_RE[i][j] = xx5;
   printf("%lf\n",C_RE[i][j]);
  }
 fclose(fC_LN);


fC_LNPN = fopen("C_LNPN.txt","r");

 for(i1 = 0; i1 < Mtc; i1++) 
   for(j1 = 0; j1 < Mtc1; j1++)
     for( i = 0; i < Mre; i++)
       for( j = 0; j < Mre1; j++){
	 fscanf(fC_LNPN,"%lf",&xx1);
	 C_RETC[i][j][i1][j1] = xx1;
//  printf("\n C_RETC=%lf",C_RETC[i][j][i1][j1]);
       }

// for(i1=0; i1<Mre; i1++)
// for(j1=0; j1<Mre1; j1++)
// for(i=0; i<Mtc; i++)
// for(j=0; j<Mtc1; j++){
//    fscanf(fC_LNPN,"%lf",&xx1);
//    C_RETC[i][j][i1][j1] = xx1;
// //    printf("\n CRETC=%lf",C_RETC[i][j][i1][j1]);
//   }
 fclose(fC_LNPN);

fC_LNLN = fopen("C_LNLN.txt","r");

for(i1 = 0; i1 < Mre; i1++)
  for(j1 = 0; j1 < Mre1; j1++)
    for(i=0; i<Mre; i++)
      for(j=0; j<Mre1; j++){
	fscanf(fC_LNLN,"%lf",&xx2);
	C_RERE[i][j][i1][j1] = xx2;
//  printf("\n C_RERE=%lf",C_RERE[i][j][i1][j1]);
      }
 fclose(fC_LNLN);


fC_PNLN = fopen("C_PNLN.txt","r");

for(i1=0; i1<Mre; i1++)
  for(j1=0; j1<Mre1; j1++)
    for(i=0; i<Mtc; i++)
      for(j=0; j<Mtc1; j++){
	fscanf(fC_PNLN,"%lf",&xx3);
	C_TCRE[i][j][i1][j1] = xx3;
// printf("\n C_TCRE=%lf",C_TCRE[i][j][i1][j1]);
      }
 fclose(fC_PNLN);

fC_PNPN = fopen("C_PNPN.txt","r");

for(i1=0; i1<Mtc; i1++)
  for(j1=0; j1<Mtc1; j1++)
    for(i=0; i<Mtc; i++)
      for(j=0; j<Mtc1; j++){
	fscanf(fC_PNPN,"%lf",&xx4);
	C_TCTC[i][j][i1][j1] = xx4;
 	printf("\n C_TCTC=%lf",C_TCTC[i][j][i1][j1]);
      }
 fclose(fC_PNPN);

/*
for(i=0; i<Mtc; i++)
for(j=0; j<Mtc1; j++){
  R = 1.0 * rand()/(RAND_MAX + 1.0);
  if(R < 0.3) {C_TC[i][j]=1; 
                printf("TC: %d %d \n", i,j);}
  }

for(i=0; i<Mre; i++)
for(j=0; j<Mre1; j++){
  R = 1.0 * rand()/(RAND_MAX + 1.0);
  if(R < 0.3) {C_RE[i][j]=1; 
                printf("RE: %d %d \n", i,j);}
  }

  
for(i=0; i<Mre; i++)
for(j=0; j<Mre1; j++)
for(i1=0; i1<Mre; i1++)
for(j1=0; j1<Mre1; j1++){
  R = 1.0 * rand()/(RAND_MAX + 1.0);
  if(R < 0.5) {C_RERE[i][j][i1][j1]=1;
//               C_RERE[i1][j1][i][j]=1;
}
  }

for(i=0; i<Mre; i++)
for(j=0; j<Mre1; j++)
for(i1=0; i1<Mtc; i1++)
for(j1=0; j1<Mtc1; j1++){
  R = 1.0 * rand()/(RAND_MAX + 1.0);
  if(R < 0.5) {C_RETC[i][j][i1][j1]=1;
//               C_TCRE[i1][j1][i][j]=1;
}
  }

for(i=0; i<Mtc; i++)
for(j=0; j<Mtc1; j++)
for(i1=0; i1<Mre; i1++)
for(j1=0; j1<Mre1; j1++){
  R = 1.0 * rand()/(RAND_MAX + 1.0);
  if(R < 0.5) {C_TCRE[i][j][i1][j1]=1;
//               C_RETC[i1][j1][i][j]=1;
}
  }

for(i=0; i<Mtc; i++)
for(j=0; j<Mtc1; j++)
for(i1=0; i1<Mtc; i1++)
for(j1=0; j1<Mtc1; j1++){
  R = 1.0 * rand()/(RAND_MAX + 1.0);
  if(R < 0.5) {C_TCTC[i][j][i1][j1]=1;
//               C_TCTC[i1][j1][i][j]=1;
}
  }
*/
/*
C_RETC[0][0][0][0]=1;
C_RETC[0][0][1][0]=1;
C_RETC[0][0][2][0]=1;
C_RETC[1][0][3][0]=1;
C_RETC[1][0][4][0]=1;
C_RETC[1][0][5][0]=1;

C_RERE[0][0][1][0]=1;
C_RERE[1][0][0][0]=1; 

C_RE[0][0]=1;
C_RE[1][0]=1;
C_TC[1][0]=1;
C_TC[3][0]=1;
*/

//--- Change Input ----------------------------------
/* srand(12);
 for(i=0; i<Mre; i++)
 for(j=0; j<Mre1; j++) {
    printf("\n %d %d %d", i,j,C_RE[i][j]);
    C_RE[i][j]=0;
 }
 for(i=0; i<Mre; i++)
 for(j=0; j<Mre1; j++){
   R = 1.0 * rand()/(RAND_MAX + 1.0);
   if(R < 0.3) {C_RE[i][j]=1;
                 printf("RE: %d %d \n", i,j);}
   }
*/

printf("\n End connections");
j=Mtc1/2;
for(i=0; i<Mtc; i++)
  printf("\n TC %ld %ld %ld", i, j, C_TC[i][j]);
j=Mtc1/4;
for(i=0; i<Mtc; i++)
  printf("\n TC %ld %ld %ld", i, j, C_TC[i][j]);
j=3*Mtc1/4;
for(i=0; i<Mtc; i++)
  printf("\n TC %ld %ld %ld", i, j, C_TC[i][j]);

for(i=0; i<Mre; i++)
  printf("\n RE %ld %ld", i, C_RE[i][Mre1/2]);

for(i1=0; i1<Mtc; i1++){
  for(i=0; i<Mre; i++) printf("%d ", C_RETC[i][0][i1][0]);
  printf("\n");
}

//------------------stimulus initialization----------------------
aP.init();
aPRE[0][0].Beta=0.0025; //0.0015;
srand(2);
   for(i=0; i < Mre; ++i)
    for(j=0; j < Mre1; ++j){
     R = 2.0 * rand()/(RAND_MAX + 1.0) - 1.0;
     if((R < -1) || (R > 1)) cout << "\n CHECK RANDOM !!!!!!!!";
     while((r[0]*r[1]>0) && (r[1]*R>0)) {R = 2.0 * rand()/(RAND_MAX + 1.0) - 1.0;}
     r[1]=r[0]; r[0]=R;
     aPRE[i][j].Beta = aPRE[0][0].Beta + R * 0.0001; 
     aPRE[0][0].Beta = 0.0025; //0.0015;
     aPRE[i][j].init();
}

//--------------open ALL files-------------------------------------
  f2 = fopen("dat", "w");
  f22 = fopen("dat1", "w");
  f222 = fopen("dat1all", "w");
  f6 = fopen("graf_re", "w");
  f7 = fopen("time_re", "w");
  f8 = fopen("graf_tc", "w");
  f9 = fopen("time_tc", "w");
  f10 = fopen("graf_cx", "w");
  f11 = fopen("time_cx", "w");
  f12 = fopen("graf_in", "w");
  f13 = fopen("time_in", "w");
  fNoise = fopen("time_Hz100Syn200", "r");
  fre = fopen("LNvariations","w");
  ftc = fopen("PNvariations","w");

//---------------show initial values of ALL variables----------------
  for(j = 0; j < Mtc1; ++j)
  for(i = 0; i < Mtc; ++i){
     if(I_TC == 1){
        printf("\n TC: ");
        for(k = 0; k < N_TC; ++k) 
           printf("%lf ", y_ini[no_tc[i][j][k]]);}
  }
  for(j = 0; j < Mre1; ++j)
  for(i = 0; i < Mre; ++i){
     if(I_RE == 1){
        printf("\n RE: ");
        for(k = 0; k < N_RE; ++k)   
           printf("%lf ", y_ini[no_re[i][j][k]]);}
  }
  for(j = 0; j < Mc1; ++j)
  for(i = 0; i < Mc; ++i){
     if(I_CX == 1){
        printf("\n CX: ");
        for(k = 0; k < N_CX; ++k) 
           printf("%lf ", y_ini[no_cx[i][j][k]]);}
     if(I_IN == 1){
        printf("\n IN: ");
        for(k = 0; k < N_IN; ++k)   
           printf("%lf ", y_ini[no_in[i][j][k]]);}
  }

  fSYNre = fopen("SYNre", "w"); 
  fSYNtc = fopen("SYNtc", "w");
  
  for(i1=0; i1<Mtc; i1++){
    fprintf(fSYNre,"%d ", i1);
    for(i=0; i<Mre; i++)
      fprintf(fSYNre,"%d ", C_RETC[i][0][i1][0]);
    fprintf(fSYNre,"\n");
  }
  for(i1=0; i1<Mtc; i1++){
    fprintf(fSYNtc,"%d ", i1);
    for(i=0; i<Mtc; i++)
      fprintf(fSYNtc,"%d ", C_TCTC[i][0][i1][0]);
    fprintf(fSYNtc,"\n");
  }
  
  fclose(fSYNre);
  fclose(fSYNtc);

//------------------------------------

 for(i=0; i<Mtc; i++)
   for(j=0; j<Mtc1; j++) 
     C_TC1[i][j]=C_TC[i][j];
 
 for(i=0; i<Mre; i++)
   for(j=0; j<Mre1; j++) 
     C_RE1[i][j]=C_RE[i][j];

tc=100000000;

//----------------CALCULATION----------------------------------------
  printf("\n CALCULATION IN PROGRESS!!!: t= %lf: tmax= %lf", t,tmax);
  ih = (int) (1/h);
  p=0;

  while( t < tmax){ 
  ii = ii + 1;       

  rk(N_EQ, fun, h, t, y_ini, f_ini, y1, y2);
  t = t + h;

if((t > 49.0)&&(t > 50.0)){
//--------the scaling of the conductances-----------------------------------
  //  g_ext3 = g_Extern_ampa/100;
  //  g_ext2 = g_Extern_ampa/100;
  g_ext1 = g_Extern_ampa1;
  g_ext2 = g_Extern_ampa2; 
  g_ext3 = g_Extern_ampa3;
  g_ext4 = g_Extern_ampa4;

if(I_RE == 1){
  g_gaba_a = g_GABA_A / re_cell[0][0].S; 
  g_ampa = g_AMPA / re_cell[0][0].S;
  g_gaba_b_re_re = g_GABA_B_RE_RE / re_cell[0][0].S;
//  g_ampa_cx_re = g_AMPA_CX_RE / re_cell[0][0].S;
  g_ext_re = g_Extern_ampa/5; ///re_cell[0][0].S;
}
if(I_TC == 1){
  g_gaba_a1 = g_GABA_A1 / tc_cell[0][0].S;
  g_ampa1 = g_AMPA1 / tc_cell[0][0].S;
  g_gaba_b = g_GABA_B / tc_cell[0][0].S;
//  g_ampa_cx_tc = g_AMPA_CX_TC / tc_cell[0][0].S;
  g_ext_tc = g_Extern_ampa; ///tc_cell[0][0].S;
}

if(I_CX == 1){
  g_ampa_cx_cx = g_AMPA_CX_CX / cx_cell[0][0].S_CX_DEND;
  g_ampa_tc_cx = g_AMPA_TC_CX / cx_cell[0][0].S_CX_DEND;
  g_gaba_a_in_cx = g_GABA_A_IN_CX / cx_cell[0][0].S_CX_DEND;
  g_gaba_b_in_cx = g_GABA_B_IN_CX / cx_cell[0][0].S_CX_DEND;
  g_ext_cx = g_Extern_ampa/cx_cell[0][0].S_CX_DEND/10; //2;
}
if(I_IN == 1){
  g_ampa_cx_in = g_AMPA_CX_IN / in_cell[0][0].S_CX_DEND;  //S_IN;
  g_ampa_tc_in = g_AMPA_TC_IN / in_cell[0][0].S_CX_DEND;  //S_IN;
//  g_ext_in = g_Extern_ampa/in_cell[0][0].S_CX_DEND / 3; // S_IN /15; //3;
}
}

Noise=Noi[ii];
//--- Change Input ----------------------------------

if((t > (tc-h/2.)) && (t < (tc+h/2.))){
 tc=tc+5000000;
 ic=ic+1;
 srand(ic);

 
 for(i=0; i<Mtc; i++)
   for(j=0; j<Mtc1; j++) {
     C_TC[i][j]=C_TC1[i][j];
 }

 for(i=0; i<Mre; i++)
   for(j=0; j<Mre1; j++) {
     C_RE[i][j]=C_RE1[i][j];
 }
/*
 fprintf(ftc,"t=%lf \n", t);
 for(i=0; i<Mtc; i++)
 for(j=0; j<Mtc1; j++){
   R = 1.0 * rand()/(RAND_MAX + 1.0);
   if((C_TC1[i][j]==0) && (R < 0.04)){
      C_TC[i][j]=1;
      fprintf(ftc,"TC+ %d %d \n", i,j);
   }}
 for(i=0; i<Mtc; i++)
 for(j=0; j<Mtc1; j++){
   R = 1.0 * rand()/(RAND_MAX + 1.0);
   if((C_TC1[i][j]==1) && (R < 0.1)){
      C_TC[i][j]=0;
      fprintf(ftc,"TC- %d %d \n", i,j);
   }}
*/
 fprintf(fre,"t=%lf \n", t);
 for(i=0; i<Mre; i++)
 for(j=0; j<Mre1; j++){
   R = 1.0 * rand()/(RAND_MAX + 1.0);
   if((C_RE1[i][j]==0) && (R < 0.2)){
      C_RE[i][j]=1;
      fprintf(fre,"RE+ %d %d \n", i);
   }}
/*
 for(i=0; i<Mre; i++)
 for(j=0; j<Mre1; j++){
   R = 1.0 * rand()/(RAND_MAX + 1.0);
   if((C_RE1[i][j]==1) && (R < 0.1)){
      C_RE[i][j]=0;
      fprintf(fre,"RE- %d %d \n", i);
   }}
*/
}


//-------------------------------------------------
if( (t>(1200-0.5*h)) && (t<(1200+0.5*h)) ){
KMAXrere=0;
KMAXtcre=0;
KMAXretc=0;
KMAXtctc=0;
KMAXretcB=0;
for(j = 0; j < Mre1; ++j)
  for(i = 0; i < Mre; ++i){
    for(k = 0; k < k_re_re[i][j]; ++k){
       KMAXrere=KMAXrere+1;
       SynapsesRERE[KMAXrere-1][0]=g_re_re[i][j][k].FF*re_cell[0][0].S;    //SynapsesRERE[KMAXrere-1][0]=g_re_re[i][j][k].FF;
    }
    for(k = 0; k < k_tc_re[i][j]; ++k){
       KMAXtcre=KMAXtcre+1;
       SynapsesTCRE[KMAXtcre-1][0]=a_tc_re[i][j][k].FF;
    }
  }
for(j = 0; j < Mtc1; ++j)
  for(i = 0; i < Mtc; ++i){
    for(k = 0; k < k_re_tc[i][j]; ++k){
       KMAXretc=KMAXretc+1;
       SynapsesRETC[KMAXretc-1][0]=g_re_tc[i][j][k].FF*tc_cell[0][0].S;         //SynapsesRETC[KMAXretc-1][0]=g_re_tc[i][j][k].FF;
    }
    for(k = 0; k < k_re_tc[i][j]; ++k){
       KMAXretcB=KMAXretcB+1;
       SynapsesRETCB[KMAXretcB-1][0]=gb[i][j][k].FF;
    }
    for(k = 0; k < k_tc_tc[i][j]; ++k){
       KMAXtctc=KMAXtctc+1;
       SynapsesTCTC[KMAXtctc-1][0]=a_tc_tc[i][j][k].FF;
    }
  }
}

if( (t>(3200-0.5*h)) && (t<(3200+0.5*h)) ){
KMAXrere=0;
KMAXtcre=0;
KMAXretc=0;
KMAXtctc=0;
KMAXretcB=0;
for(j = 0; j < Mre1; ++j)
  for(i = 0; i < Mre; ++i){
    for(k = 0; k < k_re_re[i][j]; ++k){
       KMAXrere=KMAXrere+1;
       SynapsesRERE[KMAXrere-1][1]=g_re_re[i][j][k].FF*re_cell[0][0].S;      //SynapsesRERE[KMAXrere-1][1]=g_re_re[i][j][k].FF;
    }
    for(k = 0; k < k_tc_re[i][j]; ++k){
       KMAXtcre=KMAXtcre+1;
       SynapsesTCRE[KMAXtcre-1][1]=a_tc_re[i][j][k].FF;
    }
  }
for(j = 0; j < Mtc1; ++j)
  for(i = 0; i < Mtc; ++i){
    for(k = 0; k < k_re_tc[i][j]; ++k){
       KMAXretc=KMAXretc+1;
       SynapsesRETC[KMAXretc-1][1]=g_re_tc[i][j][k].FF*tc_cell[0][0].S;        //SynapsesRETC[KMAXretc-1][1]=g_re_tc[i][j][k].FF;
    }
    for(k = 0; k < k_re_tc[i][j]; ++k){
       KMAXretcB=KMAXretcB+1;
       SynapsesRETCB[KMAXretcB-1][1]=gb[i][j][k].FF;
    }
    for(k = 0; k < k_tc_tc[i][j]; ++k){
       KMAXtctc=KMAXtctc+1;
       SynapsesTCTC[KMAXtctc-1][1]=a_tc_tc[i][j][k].FF;
    }
  }
}

if( (t>(5200-0.5*h)) && (t<(5200+0.5*h)) ){
KMAXrere=0;
KMAXtcre=0;
KMAXretc=0;
KMAXtctc=0;
KMAXretcB=0;
for(j = 0; j < Mre1; ++j)
  for(i = 0; i < Mre; ++i){
    for(k = 0; k < k_re_re[i][j]; ++k){
       KMAXrere=KMAXrere+1;
       SynapsesRERE[KMAXrere-1][2]=g_re_re[i][j][k].FF*re_cell[0][0].S;      //SynapsesRERE[KMAXrere-1][1]=g_re_re[i][j][k].FF;
    }
    for(k = 0; k < k_tc_re[i][j]; ++k){
       KMAXtcre=KMAXtcre+1;
       SynapsesTCRE[KMAXtcre-1][2]=a_tc_re[i][j][k].FF;
    }
  }
for(j = 0; j < Mtc1; ++j)
  for(i = 0; i < Mtc; ++i){
    for(k = 0; k < k_re_tc[i][j]; ++k){
       KMAXretc=KMAXretc+1;
       SynapsesRETC[KMAXretc-1][2]=g_re_tc[i][j][k].FF*tc_cell[0][0].S;;       //SynapsesRETC[KMAXretc-1][1]=g_re_tc[i][j][k].FF;
    }
    for(k = 0; k < k_re_tc[i][j]; ++k){
       KMAXretcB=KMAXretcB+1;
       SynapsesRETCB[KMAXretcB-1][2]=gb[i][j][k].FF;
    }
    for(k = 0; k < k_tc_tc[i][j]; ++k){
       KMAXtctc=KMAXtctc+1;
       SynapsesTCTC[KMAXtctc-1][2]=a_tc_tc[i][j][k].FF;
    }
  }
}

if( (t>(9200-0.5*h)) && (t<(9200+0.5*h)) ){
KMAXrere=0;
KMAXtcre=0;
KMAXretc=0;
KMAXtctc=0;
KMAXretcB=0;
for(j = 0; j < Mre1; ++j)
  for(i = 0; i < Mre; ++i){
    for(k = 0; k < k_re_re[i][j]; ++k){
       KMAXrere=KMAXrere+1;
       SynapsesRERE[KMAXrere-1][3]=g_re_re[i][j][k].FF*re_cell[0][0].S;      //SynapsesRERE[KMAXrere-1][1]=g_re_re[i][j][k].FF;
    }
    for(k = 0; k < k_tc_re[i][j]; ++k){
       KMAXtcre=KMAXtcre+1;
       SynapsesTCRE[KMAXtcre-1][3]=a_tc_re[i][j][k].FF;
    }
  }
for(j = 0; j < Mtc1; ++j)
  for(i = 0; i < Mtc; ++i){
    for(k = 0; k < k_re_tc[i][j]; ++k){
       KMAXretc=KMAXretc+1;
       SynapsesRETC[KMAXretc-1][3]=g_re_tc[i][j][k].FF*tc_cell[0][0].S;;       //SynapsesRETC[KMAXretc-1][1]=g_re_tc[i][j][k].FF;
    }
    for(k = 0; k < k_re_tc[i][j]; ++k){
       KMAXretcB=KMAXretcB+1;
       SynapsesRETCB[KMAXretcB-1][3]=gb[i][j][k].FF;
    }
    for(k = 0; k < k_tc_tc[i][j]; ++k){
       KMAXtctc=KMAXtctc+1;
       SynapsesTCTC[KMAXtctc-1][3]=a_tc_tc[i][j][k].FF;
    }
  }
}

if( (t>(15200-0.5*h)) && (t<(15200+0.5*h)) ){
KMAXrere=0;
KMAXtcre=0;
KMAXretc=0;
KMAXtctc=0;
KMAXretcB=0;
for(j = 0; j < Mre1; ++j)
  for(i = 0; i < Mre; ++i){
    for(k = 0; k < k_re_re[i][j]; ++k){
       KMAXrere=KMAXrere+1;
       SynapsesRERE[KMAXrere-1][4]=g_re_re[i][j][k].FF*re_cell[0][0].S;      //SynapsesRERE[KMAXrere-1][1]=g_re_re[i][j][k].FF;;
    }
    for(k = 0; k < k_tc_re[i][j]; ++k){
       KMAXtcre=KMAXtcre+1;
       SynapsesTCRE[KMAXtcre-1][4]=a_tc_re[i][j][k].FF;
    }
  }
for(j = 0; j < Mtc1; ++j)
  for(i = 0; i < Mtc; ++i){
    for(k = 0; k < k_re_tc[i][j]; ++k){
       KMAXretc=KMAXretc+1;
       SynapsesRETC[KMAXretc-1][4]=g_re_tc[i][j][k].FF*tc_cell[0][0].S;;       //SynapsesRETC[KMAXretc-1][1]=g_re_tc[i][j][k].FF;;
    }
    for(k = 0; k < k_re_tc[i][j]; ++k){
       KMAXretcB=KMAXretcB+1;
       SynapsesRETCB[KMAXretcB-1][4]=gb[i][j][k].FF;
    }
    for(k = 0; k < k_tc_tc[i][j]; ++k){
       KMAXtctc=KMAXtctc+1;
       SynapsesTCTC[KMAXtctc-1][4]=a_tc_tc[i][j][k].FF;
    }
  }
}

//-----------------------------------
     avr=0;
     av=0;
     for(j = 0; j < Mre1; ++j)
       for(i = 0; i < Mre; ++i)
         avr=avr + y_ini[no_re[i][j][0]];
     for(j = 0; j < Mtc1; ++j)
       for(i = 0; i < Mtc; ++i)
         av=av + y_ini[no_tc[i][j][0]];
     av = av/(Mtc1*Mtc);
     avr = avr/(Mre1*Mre);

   if((ii/(ih*1000))*(ih*1000) == ii) {
      printf("\n T= %lf ",t);
//        for(j = 0; j < Mtc1; ++j)
        for(i = 0; i < Mtc; ++i){
        if(I_TC == 1){
          printf("\n TC: ");
          for(k = 0; k < N_TC; ++k) 
             printf("%lf ", y_ini[no_tc[i][0][k]]);}
}
//        for(j = 0; j < Mre1; ++j)
        for(i = 0; i < Mre; ++i){
        if(I_RE == 1){
          printf("\n RE: ");
          for(k = 0; k < N_RE; ++k)   
             printf("%lf ", y_ini[no_re[i][0][k]]);}
}
//        for(j = 0; j < Mc1; ++j)
        for(i = 0; i < Mc; ++i){
        if(I_CX == 1){
          printf("\n CX: ");
          for(k = 0; k < N_CX; ++k) 
             printf("%lf ", cx_cell[i][0].v_SOMA);}
        if(I_IN == 1){
          printf("\n IN: ");
          for(k = 0; k < N_IN; ++k)   
             printf("%lf ", y_ini[no_in[i][j][k]]);}
}
}
   /*  
   if((t > t3D) && ((ii/(ih))*(ih) == ii)){
     for(j = 0; j < Mre1; ++j){
       for(i = 0; i < Mre; ++i){
             if(I_RE == 1) fprintf(f6,"%lf ", y_ini[no_re[i][j][0]]);}
       fprintf(f6,"\n");}

     for(j = 0; j < Mtc1; ++j){
       for(i = 0; i < Mtc; ++i){
             if(I_TC == 1) fprintf(f8,"%lf ", y_ini[no_tc[i][j][0]]); }
       fprintf(f8,"\n");}

     for(j = 0; j < Mc1; ++j){
       for(i = 0; i < Mc; ++i){
//         fprintf(f10,"%lf ", y_ini[no_cx[i][j][0]]); 
         if(I_CX == 1) fprintf(f10,"%lf ", cx_cell[i][j].v_SOMA); 
         if(I_IN == 1) fprintf(f12,"%lf ", in_cell[i][j].v_SOMA); 
         }
       fprintf(f10,"\n");
       fprintf(f12,"\n");
       }
  }
  */
   if((t > ttime) && ((ii/25)*(25) == ii)){
     if(I_RE == 1) fprintf(f7,"%lf %lf %lf ", t, avr,
                             aPRE[0][0].RS*a_ext1.g/re_cell[0][0].S);
			   //                  gb_re_re[0][0][0].G);
     if(I_TC == 1) fprintf(f9,"%lf %lf %lf %lf %lf ", t, av,
                             aP.RS*a_ext1.g/tc_cell[0][0].S,
                             gb[0][0][0].G, g_re_tc[0][0][0].G);
     for(i = 0; i < Mre; ++i){
         if(I_RE == 1) fprintf(f7,"%lf ", y_ini[no_re[i][0][0]]);}
     fprintf(f7,"\n");

     for(i = 0; i < Mtc; ++i){
         if(I_TC == 1) fprintf(f9,"%lf ", y_ini[no_tc[i][0][0]]);}
     fprintf(f9,"\n");
   }

   if((t > ttime) && ((ii/(2))*(2) == ii)){
     if(I_CX == 1) fprintf(f11,"%lf ", t);
     if(I_IN == 1) fprintf(f13,"%lf ", t);
     for(i = 0; i < Mc; ++i){
         if(I_CX == 1) fprintf(f11,"%lf ", cx_cell[i][Mc1/2].v_SOMA);
         if(I_IN == 1) fprintf(f13,"%lf ", in_cell[i][Mc1/2].v_SOMA);
     }
     fprintf(f11,"\n");
     fprintf(f13,"\n");
   }

}

fsRERE = fopen("SYnapceRERE", "w");
fsRETC = fopen("SYnapceRETC", "w");
fsTCRE = fopen("SYnapceTCRE", "w");
fsTCTC = fopen("SYnapceTCTC", "w");
fsRETCB = fopen("SYnapceRETCB", "w"); 

for(k=0; k < KMAXrere; k++){
   for(j=0; j < 5; j++)
      fprintf(fsRERE,"%lf ",SynapsesRERE[k][j]);
   fprintf(fsRERE,"\n ");
} 
for(k=0; k < KMAXretc; k++){
   for(j=0; j < 5; j++)
      fprintf(fsRETC,"%lf ",SynapsesRETC[k][j]);
   fprintf(fsRETC,"\n ");
}
for(k=0; k < KMAXretcB; k++){
   for(j=0; j < 5; j++)
      fprintf(fsRETCB,"%lf ",SynapsesRETCB[k][j]);
   fprintf(fsRETCB,"\n ");
}
for(k=0; k < KMAXtcre; k++){
   for(j=0; j < 5; j++)
      fprintf(fsTCRE,"%lf ",SynapsesTCRE[k][j]);
   fprintf(fsTCRE,"\n ");
}
for(k=0; k < KMAXtctc; k++){
   for(j=0; j < 5; j++)
      fprintf(fsTCTC,"%lf ",SynapsesTCTC[k][j]);
   fprintf(fsTCTC,"\n ");
}

fclose(fsRERE);
fclose(fsRETC);
fclose(fsTCRE);
fclose(fsTCTC);
fclose(fsRETCB);
//--------------------END CALCULATION-------------------------------

//-----------------close ALL files-----------------------------------
  fclose(f2);
  fclose(f22);
  fclose(f6);
  fclose(f7);
  fclose(f8);
  fclose(f9);
  fclose(f10);
  fclose(f11);
  fclose(f12);
  fclose(f13);
  //  fclose(fNoise);
  fclose(fre); 
  printf("\n"); 
}

//+++++++++++ Function to calculate the index of BOUNDARY elements ++++++++++
int b(unsigned btype, unsigned m, int ind) {
      if(btype == 0) {  //------flow boundary conditions:
        if( (ind >= 0) && (ind <= (m-1)) ) return( ind );
        if(ind < 0) return( -ind );
        if(ind > (m-1)) return( 2*m - ind - 2); }

      else if(btype == 1) {  //------periodic boundary conditions: 
        if( (ind >= 0) && (ind <= (m-1)) ) return( ind );
        if(ind < 0) return( ind + m);
        if(ind > (m-1)) return (ind - m); }

      else if(btype == 9) {  //------NO boundary conditions:
        return( ind ); }     
}

//+++++++++++ Function to calculate the right sides for ALL ODE +++++++++++++++
void fun(unsigned neq, double x, double *y_ini, double *f_ini){
double scale;
int i, j, k, i1, j1, ii, jj, nst, kk[Max][Max1];

//========here the MAIN loop to calculate intrinsic conductances===========
//--------(f_ini IS changed, y_ini IS NOT changed)-------------------------
for(i=0; i < Mre; ++i)
for(j=0; j < Mre1; ++j){
  if(I_RE == 1) re_cell[i][j].calc(x, y_ini+no_re[i][j][0], 
                                      f_ini+no_re[i][j][0]); }

for(i=0; i < Mtc; ++i)
for(j=0; j < Mtc1; ++j){
  if(I_TC == 1) tc_cell[i][j].calc(x, y_ini+no_tc[i][j][0], 
                                      f_ini+no_tc[i][j][0]); }

for(i=0; i < Mc; ++i)
for(j=0; j < Mc1; ++j){
  if(I_CX == 1) cx_cell[i][j].calc(x, y_ini+no_cx[i][j][0], 
                                      f_ini+no_cx[i][j][0]);
  if(I_IN == 1) in_cell[i][j].calc(x, y_ini+no_in[i][j][0], 
                                      f_ini+no_in[i][j][0]); }

//========here the MAIN loop to calculate synaptic conductances=============
//--------(f_ini IS changed, y_ini IS NOT changed) -------------------------
for(i = 0; i < Mre; ++i)
for(j = 0; j < Mre1; ++j){

//--------reciprocal GABA-A between RE cells---------------------------------
if(I_RE == 1){
  for(i1 = 0, k = 0; i1 < Mre; ++i1)
  for(j1 = 0; j1 < Mre1; ++j1){

//    scale = sqrt((double) (i1*i1 + j1*j1));
   if(C_RERE[i1][j1][i][j] > 0){
    ii = i1; jj = j1;
//  if( (BOUND != 9) || ( (ii>=0) && (ii<=(Mre-1)) && 
//                                   (jj>=0) && (jj<=(Mre1-1)) ) ){
          if( !( (ii == i) && (jj == j) && (SELFING == 0) ) ){
              g_re_re[i][j][k].calc(C_RERE[i1][j1][i][j]/re_cell[0][0].S, x, y_ini+no_ga_rere[i][j][k][0],    //g_re_re[i][j][k].calc(g_gaba_a*C_RERE[i1][j1][i][j], x, y_ini+no_ga_rere[i][j][k][0],
                      f_ini+no_ga_rere[i][j][k][0],
                      y_ini[no_re[i][j][0]], y_ini[no_re[ii][jj][0]]);
          if(I_GB_RERE == 1){ gb_re_re[i][j][k].calc(g_gaba_b_re_re, x, 
                   y_ini+no_gb_rere[i][j][k][0], f_ini+no_gb_rere[i][j][k][0], 
                   y_ini[no_re[i][j][0]], y_ini[no_re[ii][jj][0]]); }
              ++k; }   // number of synapces of the cell i-j 
// }
  } }
  k_re_re[i][j] = k;
  for(k = 0; k < k_re_re[i][j]; ++k){
    if(k_re_re[i][j] == 0) { }  //NO synapces
    else { 
       g_re_re[i][j][k].I = g_re_re[i][j][k].I / k_re_re[i][j];
       f_ini[no_re[i][j][0]] = f_ini[no_re[i][j][0]] - g_re_re[i][j][k].I; 
      if(I_GB_RERE == 1) {gb_re_re[i][j][k].I = gb_re_re[i][j][k].I / 
                                                k_re_re[i][j];
                          f_ini[no_re[i][j][0]] = f_ini[no_re[i][j][0]] - 
			  gb_re_re[i][j][k].I; }
}
}}

//-----------AMPA from TC to RE cells---------------------------------------
if(I_RE == 1 && I_TC == 1){
  for(i1 = 0, k = 0; i1 < Mtc; ++i1)
  for(j1 = 0; j1 < Mtc1; ++j1){
//    scale = sqrt((double) (i1*i1 + j1*j1));
//    if(scale <= MS_TC_RE_MAX){
   if(C_TCRE[i1][j1][i][j] > 0){
    ii = i1; jj = j1;
//  if( (BOUND != 9) || ( (ii>=0) && (ii<=(Mtc-1)) && 
//                                   (jj>=0) && (jj<=(Mtc1-1)) ) ){
        a_tc_re[i][j][k].calc(g_ampa, x, y_ini[no_re[i][j][0]], 
                           y_ini[no_tc[ii][jj][0]]);
        ++k;  // number of synapces of the cell i-j 
// }
  } }
  k_tc_re[i][j] = k;
  for(k = 0; k < k_tc_re[i][j]; ++k){
    if(k_tc_re[i][j] == 0) { }  //NO synapces
    else {
       a_tc_re[i][j][k].I = a_tc_re[i][j][k].I / k_tc_re[i][j];
       f_ini[no_re[i][j][0]] = f_ini[no_re[i][j][0]] - a_tc_re[i][j][k].I; }
}}

//-------------AMPA from CX to RE cells-------------------------------------
if(I_CX == 1 && I_RE == 1){
  for(i1 = -MS_CX_TC, k = 0; i1 <= MS_CX_TC; ++i1)
  for(j1 = -MS_CX_TC1; j1 <= MS_CX_TC1; ++j1){
    scale = sqrt((double) (i1*i1 + j1*j1));
    if(scale <= MS_CX_TC_MAX){
      ii = i + i1; jj = j + j1;
      if( (BOUND != 9) || ( (ii>=0) && (ii<=(Mc-1)) && 
                                       (jj>=0) && (jj<=(Mc1-1)) ) ){
        a_cx_re[i][j][k].calc(g_ampa_cx_re, x, y_ini[no_re[i][j][0]], 
                           cx_cell[b(BOUND,Mc,ii)][b(BOUND,Mc1,jj)].v_SOMA);
        ++k; } // number of synapces of the cell i-j 
  } }
  kk[i][j] = k;
  for(k = 0; k < kk[i][j]; ++k){
    if(kk[i][j] == 0) { }  //NO synapces
    else {
       a_cx_re[i][j][k].I = a_cx_re[i][j][k].I / kk[i][j];
       f_ini[no_re[i][j][0]] = f_ini[no_re[i][j][0]] - a_cx_re[i][j][k].I; }
}}
}


for(i = 0; i < Mtc; ++i)
for(j = 0; j < Mtc1; ++j){

//--------reciprocal AMPA between TC cells----------------------------------
 if(I_TC == 1){
  for(i1 = 0, k = 0; i1 < Mtc; ++i1)
  for(j1 = 0; j1 < Mtc1; ++j1){

//    scale = sqrt((double) (i1*i1 + j1*j1));
//    if(scale <= MS_TC_TC_MAX){
   if(C_TCTC[i1][j1][i][j] > 0){
    ii = i1; jj = j1;
//  if( (BOUND != 9) || ( (ii>=0) && (ii<=(Mtc-1)) && 
//                                     (jj>=0) && (jj<=(Mtc1-1)) ) ){
          if( !( (ii == i) && (jj == j) && (SELFING == 0) ) ){
              a_tc_tc[i][j][k].calc(g_ampa1, x,  y_ini[no_tc[i][j][0]], 
                               y_ini[no_tc[ii][jj][0]]);
              ++k; }   // number of synapces of the cell i-j 
// }
  } }
  k_tc_tc[i][j] = k;
  for(k = 0; k < k_tc_tc[i][j]; ++k){
    if(k_tc_tc[i][j] == 0) { }  //NO synapces
    else { 
       a_tc_tc[i][j][k].I = a_tc_tc[i][j][k].I / k_tc_tc[i][j];
       f_ini[no_tc[i][j][0]] = f_ini[no_tc[i][j][0]] - a_tc_tc[i][j][k].I; }

}}

//---------GABA-A and GABA-B from RE to TC cells------------------------------
if(I_RE == 1 && I_TC == 1){
  for(i1 = 0, k = 0; i1 < Mre; ++i1)
  for(j1 = 0; j1 < Mre1; ++j1){
//    scale = sqrt((double) (i1*i1 + j1*j1));
//    if(scale <= MS_RE_TC_MAX){
   if(C_RETC[i1][j1][i][j] > 0){
    ii = i1; jj = j1;
//  if( (BOUND != 9) || ( (ii>=0) && (ii<=(Mre-1)) && 
//                                     (jj>=0) && (jj<=(Mre1-1)) ) ){
       g_re_tc[i][j][k].calc(C_RETC[i1][j1][i][j]/tc_cell[0][0].S, x, y_ini+no_ga_retc[i][j][k][0],    //g_re_tc[i][j][k].calc(g_gaba_a1*C_RETC[i1][j1][i][j], x, y_ini+no_ga_retc[i][j][k][0],
                      f_ini+no_ga_retc[i][j][k][0],
                      y_ini[no_tc[i][j][0]], y_ini[no_re[ii][jj][0]]);
       if(I_GB == 1){ gb[i][j][k].calc(g_gaba_b, x, y_ini+no_g[i][j][k][0],
                           f_ini+no_g[i][j][k][0], y_ini[no_tc[i][j][0]], 
                             y_ini[no_re[ii][jj][0]]); }
         ++k;    // number of synapces of the cell i-j
// }
  } }
  k_re_tc[i][j] = k;
  for(k = 0; k < k_re_tc[i][j]; ++k){
    if(k_re_tc[i][j] == 0) { }  //NO synapces
    else {
       g_re_tc[i][j][k].I = g_re_tc[i][j][k].I / k_re_tc[i][j];
       if(I_GB == 1) gb[i][j][k].I = gb[i][j][k].I / k_re_tc[i][j];
       f_ini[no_tc[i][j][0]] = f_ini[no_tc[i][j][0]] - g_re_tc[i][j][k].I;
       if(I_GB == 1) f_ini[no_tc[i][j][0]] = f_ini[no_tc[i][j][0]] - 
                                                          gb[i][j][k].I; 
}
}}

//------------AMPA from CX to TC cells-------------------------------------
if(I_CX == 1 && I_TC == 1){
  for(i1 = -MS_CX_TC, k = 0; i1 <= MS_CX_TC; ++i1)
  for(j1 = -MS_CX_TC1; j1 <= MS_CX_TC1; ++j1){
    scale = sqrt((double) (i1*i1 + j1*j1));
    if(scale <= MS_CX_TC_MAX){
      ii = i + i1; jj = j + j1;
      if( (BOUND != 9) || ( (ii>=0) && (ii<=(Mc-1)) && 
                                       (jj>=0) && (jj<=(Mc1-1)) ) ){
        a_cx_tc[i][j][k].calc(g_ampa_cx_tc, x, y_ini[no_tc[i][j][0]], 
                           cx_cell[b(BOUND,Mc,ii)][b(BOUND,Mc1,jj)].v_SOMA);
        ++k; } // number of synapces of the cell i-j 
  } }
  kk[i][j] = k;
  for(k = 0; k < kk[i][j]; ++k){
    if(kk[i][j] == 0) { }  //NO synapces
    else {
       a_cx_tc[i][j][k].I = a_cx_tc[i][j][k].I / kk[i][j];
       f_ini[no_tc[i][j][0]] = f_ini[no_tc[i][j][0]] - a_cx_tc[i][j][k].I; }
}}
}

for(i = 0; i < Mc; ++i)
for(j = 0; j < Mc1; ++j){

//-----------reciprocal AMPA between CX cells--------------------------------
if(I_CX == 1){
  for(i1 = -MS_CX_CX, k = 0; i1 <= MS_CX_CX; ++i1)
  for(j1 = -MS_CX_CX1; j1 <= MS_CX_CX1; ++j1){
    scale = sqrt((double) (i1*i1 + j1*j1));
    if(scale <= MS_CX_CX_MAX){
      ii = i + i1; jj = j + j1;
      if( (BOUND != 9) || ( (ii>=0) && (ii<=(Mc-1)) && 
                                       (jj>=0) && (jj<=(Mc1-1)) ) ){
        if( !( ((i1*i1 + j1*j1) == 0) && (SELFINGcx == 0) ) ){
        a_cx_cx[i][j][k].calc(g_ampa_cx_cx, x, y_ini[no_cx[i][j][0]], 
                           cx_cell[b(BOUND,Mc,ii)][b(BOUND,Mc1,jj)].v_SOMA);
        ++k; }} // number of synapces of the cell i-j 
  } }
  kk[i][j] = k;
  for(k = 0; k < kk[i][j]; ++k){
    if(kk[i][j] == 0) { }  //NO synapces
    else {
       a_cx_cx[i][j][k].I = a_cx_cx[i][j][k].I / kk[i][j];
       f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - a_cx_cx[i][j][k].I; }
}}

//-------------AMPA from CX to IN cells--------------------------------------
if(I_CX == 1 && I_IN == 1){
  for(i1 = -MS_CX_IN, k = 0; i1 <= MS_CX_IN; ++i1)
  for(j1 = -MS_CX_IN1; j1 <= MS_CX_IN1; ++j1){
    scale = sqrt((double) (i1*i1 + j1*j1));
    if(scale <= MS_CX_IN_MAX){
      ii = i + i1; jj = j + j1;
      if( (BOUND != 9) || ( (ii>=0) && (ii<=(Mc-1)) && 
                                       (jj>=0) && (jj<=(Mc1-1)) ) ){
        a_cx_in[i][j][k].calc(g_ampa_cx_in, x, y_ini[no_in[i][j][0]], 
                           cx_cell[b(BOUND,Mc,ii)][b(BOUND,Mc1,jj)].v_SOMA);
        ++k; } // number of synapces of the cell i-j 
  } }
  kk[i][j] = k;
  for(k = 0; k < kk[i][j]; ++k){
    if(kk[i][j] == 0) { }  //NO synapces
    else {
       a_cx_in[i][j][k].I = a_cx_in[i][j][k].I / kk[i][j];
       f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] - a_cx_in[i][j][k].I; }
}}

//--------------GABA-A from IN to CX cells-----------------------------------
if(I_CX == 1 && I_IN == 1){
  for(i1 = -MS_IN_CX, k = 0; i1 <= MS_IN_CX; ++i1)
  for(j1 = -MS_IN_CX1; j1 <= MS_IN_CX1; ++j1){
    scale = sqrt((double) (i1*i1 + j1*j1));
    if(scale <= MS_IN_CX_MAX){
      ii = i + i1; jj = j + j1;
      if( (BOUND != 9) || ( (ii>=0) && (ii<=(Mc-1)) && 
                                       (jj>=0) && (jj<=(Mc1-1)) ) ){
        ga_in_cx[i][j][k].calc(g_gaba_a_in_cx, x, y_ini[no_cx[i][j][0]], 
                           in_cell[b(BOUND,Mc,ii)][b(BOUND,Mc1,jj)].v_SOMA);
//                           y_ini[no_in[b(BOUND,Mc,ii)][b(BOUND,Mc1,jj)][0]]);
        ++k; } // number of synapces of the cell i-j 
  } }
  kk[i][j] = k;
  for(k = 0; k < kk[i][j]; ++k){
    if(kk[i][j] == 0) { }  //NO synapces
    else {
       ga_in_cx[i][j][k].I = ga_in_cx[i][j][k].I / kk[i][j];
       f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - ga_in_cx[i][j][k].I; }
}}

//--------------GABA-B from IN to CX cells----------------------------------
if(I_CX == 1 && I_IN == 1){
  for(i1 = -MS_IN_CX, k = 0; i1 <= MS_IN_CX; ++i1)
  for(j1 = -MS_IN_CX1; j1 <= MS_IN_CX1; ++j1){
    scale = sqrt((double) (i1*i1 + j1*j1));
    if(scale <= MS_IN_CX_MAX){
      ii = i + i1; jj = j + j1;
      if( (BOUND != 9) || ( (ii>=0) && (ii<=(Mc-1)) && 
                                       (jj>=0) && (jj<=(Mc1-1)) ) ){
        gb_in_cx[i][j][k].calc(g_gaba_b_in_cx, x, y_ini+no_gcx[i][j][k][0],
                           f_ini+no_gcx[i][j][k][0], y_ini[no_cx[i][j][0]], 
                           in_cell[b(BOUND,Mc,ii)][b(BOUND,Mc1,jj)].v_SOMA);
//                           y_ini[no_in[b(BOUND,Mc,ii)][b(BOUND,Mc1,jj)][0]]);
        ++k; } // number of synapces of the cell i-j 
  } }
  kk[i][j] = k;
  for(k = 0; k < kk[i][j]; ++k){
    if(kk[i][j] == 0) { }  //NO synapces
    else {
       gb_in_cx[i][j][k].I = gb_in_cx[i][j][k].I / kk[i][j];
       f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - gb_in_cx[i][j][k].I; }
}}

//--------------AMPA from TC to CX cells-------------------------------------
if(I_CX == 1 && I_TC == 1){
  for(i1 = -MS_TC_CX, k = 0; i1 <= MS_TC_CX; ++i1)
  for(j1 = -MS_TC_CX1; j1 <= MS_TC_CX1; ++j1){
    scale = sqrt((double) (i1*i1 + j1*j1));
    if(scale <= MS_TC_CX_MAX){
      ii = i + i1; jj = j + j1;
      if( (BOUND != 9) || ( (ii>=0) && (ii<=(Mtc-1)) && 
                                       (jj>=0) && (jj<=(Mtc1-1)) ) ){        
        a_tc_cx[i][j][k].calc(g_ampa_tc_cx, x, y_ini[no_cx[i][j][0]], 
                           y_ini[no_tc[b(BOUND,Mtc,ii)][b(BOUND,Mtc1,jj)][0]]);
        ++k; } // number of synapces of the cell i-j 
  } }
  kk[i][j] = k;
  for(k = 0; k < kk[i][j]; ++k){
    if(kk[i][j] == 0) { }  //NO synapces
    else {
       a_tc_cx[i][j][k].I = a_tc_cx[i][j][k].I / kk[i][j];
       f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - a_tc_cx[i][j][k].I; }
}}
//-------------AMPA from TC to IN cells--------------------------------------
if(I_IN == 1 && I_TC == 1){
  for(i1 = -MS_TC_IN, k = 0; i1 <= MS_TC_IN; ++i1)
  for(j1 = -MS_TC_IN1; j1 <= MS_TC_IN1; ++j1){
    scale = sqrt((double) (i1*i1 + j1*j1));
    if(scale <= MS_TC_IN_MAX+0.0001){
      ii = i + i1; jj = j + j1;
      if( (BOUND != 9) || ( (ii>=0) && (ii<=(Mtc-1)) && 
                                       (jj>=0) && (jj<=(Mtc1-1)) ) ){
       a_tc_in[i][j][k].calc(g_ampa_tc_in, x, y_ini[no_in[i][j][0]], 
                           y_ini[no_tc[b(BOUND,Mtc,ii)][b(BOUND,Mtc1,jj)][0]]);
        ++k; } // number of synapces of the cell i-j 
  } }
  kk[i][j] = k;
  for(k = 0; k < kk[i][j]; ++k){
    if(kk[i][j] == 0) { }  //NO synapces
    else {
       a_tc_in[i][j][k].I = a_tc_in[i][j][k].I / kk[i][j];
       f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] - a_tc_in[i][j][k].I; }
}}

}
//=============ENF of MAIN loop==================================================

//------------------expernal stimulation---------------------------------
nst = 0;   //type of external stimulation

//-----------train of shocks----------------------------------------
if(nst == 0){ 

  for(i = 0; i < Mtc; ++i)
    for(j = 0; j < Mtc1; ++j){
      a_ext3.calc(g_ext3, x);
      if(I_TC == 1) f_ini[no_tc[i][j][0]] = f_ini[no_tc[i][j][0]] - 
	                     y_ini[no_tc[i][j][0]] * a_ext3.g/tc_cell[0][0].S;
    }

  for(i = 0; i < Mre; ++i)
    for(j = 0; j < Mre1; ++j){
      a_ext4.calc(g_ext4, x);
      if(I_RE == 1) f_ini[no_re[i][j][0]] = f_ini[no_re[i][j][0]] - 
	                     y_ini[no_re[i][j][0]] * a_ext4.g/re_cell[0][0].S;
    }

//      if(g_ext1 > 1.e-10){
//      fscanf(fNoise,"%lf %lf",&tt, &Noise);
  a_ext1.Noise = Noise;
  a_ext2.Noise = Noise;
      a_ext1.calc(g_ext1, x);
      a_ext2.calc(g_ext2, x);

    aP.calc(x);


      
//    if(I_TC == 1) a_ext1.calc(g_ext_tc, x);
//    if(I_RE == 1) a_ext2.calc(g_ext_re, x);
//    if(I_CX == 1) a_ext3.calc(g_ext_cx, x);
//    if(I_IN == 1) a_ext4.calc(g_ext_in, x);

      for(i = 0; i < Mc; ++i)
      for(j = 0; j < Mc1; ++j){
         if(I_IN == 1) f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] - 
                       y_ini[no_in[i][j][0]] * a_ext4.g; 
         if(I_CX == 1) f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - 
                       y_ini[no_cx[i][j][0]] * a_ext3.g; 
      }
      for(i = 0; i < Mtc; ++i)
      for(j = 0; j < Mtc1; ++j){
        if(C_TC[i][j] > 0)
           if(I_TC == 1) f_ini[no_tc[i][j][0]] = f_ini[no_tc[i][j][0]] - 
	             C_TC[i][j]*y_ini[no_tc[i][j][0]] * aP.RS * a_ext1.g/tc_cell[0][0].S;
           }

      for(i = 0; i < Mre; ++i)
      for(j = 0; j < Mre1; ++j){
         aPRE[i][j].calc(x);
         if(C_RE[i][j] > 0){
           if(I_RE == 1) f_ini[no_re[i][j][0]] = f_ini[no_re[i][j][0]] - 
	     C_RE[i][j]*y_ini[no_re[i][j][0]] * aPRE[i][j].RS * a_ext2.g/re_cell[0][0].S;}
	   } 
      //     }  
} 
















//-----------ONE shock--------------------------------------------- 
else if(nst == 1) {
  if( ((x > 2000) && (x < 2050)) || ((x > 9000) && (x < 9050)) ) {
     a_ext1.calc(g_ext_tc, x);
  
        if(I_TC == 1) f_ini[no_tc[0][0][0]] = f_ini[no_tc[0][0][0]] -
           a_ext1.g * y_ini[no_tc[0][0][0]];
  }
}

//-----------temporal stuff------------------------------------------
//if(x > 1000) cx_cell[0][0].I_Stim1 = 0.375; //0.5 would correspond to 100 for soma
//if(x > 2000) cx_cell[0][0].I_Stim1 = 0.;

//if(x > 1000) cx_cell[0][0].I_Stim2 = 70;
//if(x > 2000) cx_cell[0][0].I_Stim2 = 0;

//if(x > 1000) cx_cell[0][0].I_Stim2 = 600;
//if(x > 1020) cx_cell[0][0].I_Stim2 = 0;

//if((x > 1000) && (x < 2000)) f_ini[no_cx[0][0][0]] = f_ini[no_cx[0][0][0]] + 
//                                                                      (0.666666);
//if((x > 3000) && (x < 4000)) cx_cell[0][0].v_SOMA = cx_cell[0][0].v_SOMA + (1.);
//if((x > 1900) && (x < 2000)) re_cell[0][0].f[0] = re_cell[0][0].f[0] - (1.3);
//if((x > 1000) && (x < 1200)) y_ini[no_in[0][0][0]] = y_ini[no_in[0][0][0]] + (0.5);
//if((x > 7000) && (x < 7200)) re_cell[0][0].f[0] = re_cell[0][0].f[0] - (0.5);
//if((x > 0) && (x < 100)) tc_cell[0][0].f[0] = tc_cell[0][0].f[0] - (0.5);
//if((x > 2000) && (x < 2100)) tc_cell[0][0].f[0] = tc_cell[0][0].f[0] + (0.3);
//if((x > 1200) && (x < 1250)) tc_cell[0][0].f[0] = tc_cell[0][0].f[0] - (0.3);
//if((x > 1440) && (x < 1500)) tc_cell[0][0].f[0] = tc_cell[0][0].f[0] - (0.3);
//if((x > 1700) && (x < 1800)) tc_cell[0][0].f[0] = tc_cell[0][0].f[0] - (0.3);
//if((x > 10000) && (x < 11000)) tc_cell[0][0].f[0] = tc_cell[0][0].f[0] - (0.2);
//if((x > 15000) && (x < 15500)) tc_cell[0][0].f[0] = tc_cell[0][0].f[0] + (1.7);
//if((x > 2500) && (x < 3000)) tc_cell[0][0].f[0] = tc_cell[0][0].f[0] + (0.2);
}

/*
void solout (long nr, double xold, double x, double* y, unsigned n, int* irtrn) 
{ 
  static double xout, xout1; 
  char format99[] = "x=%f  y=%12.10f %12.10f  nstep=%li\r\n"; 
  int i, j, k;

  if (nr == 1) {  
//    printf ( "x=%f  y=%12.10f %12.10f  nstep=%li\r\n", x, y[0], y[1], nr-1); 

    fprintf(f6,"%lf ", x);
    fprintf(f8,"%lf ", x);
    for(j = 0; j < M1; ++j){
       for(i = 0; i < M; ++i){
         fprintf(f6,"%lf ", y[no_re[i][j]]);
         fprintf(f8,"%lf ", y[no_tc[i][j]]); }
    fprintf(f6,"\n");
    fprintf(f8,"\n");}
       
    xout = x + 0.2; 
  } 
  else {
    while (x >= xout1) {
      printf("\n t=%lf", xout1);
      xout1 += 100.0; } 

    while (x >= xout) { 
//      printf (format99, xout, contd8(0,xout), contd8(1,xout), nr-1); 

      fprintf(f6,"%lf ", xout);
      fprintf(f8,"%lf ", xout);
      for(j = 0; j < M1; ++j){
         for(i = 0; i < M; ++i){
           fprintf(f6,"%lf ", contd8(no_re[i][j],xout) );
           fprintf(f8,"%lf ", contd8(no_tc[i][j],xout) ); }
      fprintf(f6,"\n");
      fprintf(f8,"\n");}

      xout += 0.2; 
    } 
  }
} */













