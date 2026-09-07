// Cardiac model ---------- ToR-ORd

// New adding I_NaL in HRd Model
// times or  multiple PARAs
// CaL, Ks, Kr and NaL

/*
        Cardiac model ToR-ORd
%     Copyright (C) 2019 Jakub Tomek. Contact: jakub.tomek.mff@gmail.com
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>// TorORd model
*/

// function [ time, X, parameters ] = modelRunner( X0, options, parameters, beats, ignoreFirst)
/*% A function for simulating models with various parameters. It serves as
% the interface between user's scripts (defining parameters) and simulation
% core. Here, the structure of parameters is unpacked and passed to the model
% INPUTS:
% X0 - starting state.

% options - options to ode15s (may be left empty if not
sure what to do with it).

% parameters - passed to the model and specify,
e.g. scaling of conductances or which implementation of
a current is used.

% beats - the number of beats in the simulation

% ignoreFirst - specifies how many action potentials are to
be ignored out of beats. This serves to save memory

% OUTPUTS:
% time - a cell array; the i-th element gives the timeline of the i-th
% action potential

% X - a cell array; the i-th element gives the matrix of state variables (#rows = length of
% time vector, columns are state variables).
%
% parameters - the same parameters that were passed to this function

%% Parameters are set here
% defaults which may be overwritten*/

// Cell type (endo/epi/mid)
#define         celltype            0

// Extracellular concentrations
#define         nao                 140.0
#define         cao                 1.8
#define         ko                  5.0

//Localization of ICaL and NCX: the fraction in junctional subspace
#define         ICaL_fractionSS     0.8
#define         INaCa_fractionSS    0.35


//function output=model_Torord(t,X,flag_ode, cellType, ICaL_Multiplier, ...
    //INa_Multiplier, Ito_Multiplier, INaL_Multiplier, IKr_Multiplier, IKs_Multiplier, IK1_Multiplier, IKb_Multiplier,INaCa_Multiplier,...
    //INaK_Multiplier, INab_Multiplier, ICab_Multiplier, IpCa_Multiplier, ICaCl_Multiplier, IClb_Multiplier, Jrel_Multiplier,Jup_Multiplier, nao,cao,ko,ICaL_fractionSS,INaCa_fractionSS, stimAmp, stimDur, vcParameters, apClamp, extraParams)

// physical constants
#define         R                   8314.0
#define         T                   310.0
#define         F                   96485.0

// cell geometry
#define         L1                  (0.01)
#define         rad                 (0.0011)
#define         vcell               (1000.0*3.14*rad*rad*L1)
#define         Ageo                (2.0*3.14*rad*rad+2.0*3.14*rad*L1)
#define         Acap                (2.0*Ageo)
#define         vmyo                (0.68*vcell)
#define         vnsr                (0.0552*vcell)
#define         vjsr                (0.0048*vcell)
#define         vss                 (0.02*vcell)

#define         cli                 24.0   								// Intracellular Cl  [mM]
#define         clo                 150.0  								// Extracellular Cl  [mM]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//CaMK constants
#define         KmCaMK              0.15

#define         aCaMK               0.05
#define         bCaMK               0.00068
#define         CaMKo               0.05
#define         KmCaM               0.0015

//update CaMK
//#define       CaMKb               (CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass))
//#define       CaMKa               (CaMKb+CaMKt)

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//reversal potentials
#define         ENa                 ((R*T/F)*log(nao/nai))
#define         EK                  ((R*T/F)*log(ko/ki))
#define         PKNa                0.01833
#define         EKs                 ((R*T/F)*log((ko+PKNa*nao)/(ki+PKNa*nai)))

//convenient shorthand calculations
#define         vffrt               (v*F*F/(R*T))
#define         vfrt                (v*F/(R*T))
#define         frt                 (F/(R*T))

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define         fINap               ((1.0/(1.0+KmCaMK/CaMKa)))
#define         fINaLp              ((1.0/(1.0+KmCaMK/CaMKa)))
#define         fItop               ((1.0/(1.0+KmCaMK/CaMKa)))
#define         fICaLp              ((1.0/(1.0+KmCaMK/CaMKa)))

#define 		pos(var_th,x,y) ((Length)*NDIM*(y)+NDIM*(x)+var_th)

/* Aliases of state variables. */
#define         st                  state_var[pos(0,idx,idy)]

//give names to the state vector values
#define         v                   state_var[pos(1,idx,idy)]
#define         nai                 state_var[pos(2,idx,idy)]
#define         nass                state_var[pos(3,idx,idy)]
#define         ki                  state_var[pos(4,idx,idy)]
#define         kss                 state_var[pos(5,idx,idy)]
#define         cai                 state_var[pos(6,idx,idy)]
#define         cass                state_var[pos(7,idx,idy)]
#define         cansr               state_var[pos(8,idx,idy)]
#define         cajsr               state_var[pos(9,idx,idy)]
#define         m                   state_var[pos(10,idx,idy)]
#define         hp                  state_var[pos(11,idx,idy)]
#define         h                   state_var[pos(12,idx,idy)]
#define         j                   state_var[pos(13,idx,idy)]

#define         jp                  state_var[pos(14,idx,idy)]
#define         mL                  state_var[pos(15,idx,idy)]
#define         hL                  state_var[pos(16,idx,idy)]
#define         hLp                 state_var[pos(17,idx,idy)]
#define         a                   state_var[pos(18,idx,idy)]
#define         iF                  state_var[pos(19,idx,idy)]
#define         iS                  state_var[pos(20,idx,idy)]
#define         ap                  state_var[pos(21,idx,idy)]
#define         iFp                 state_var[pos(22,idx,idy)]
#define         iSp                 state_var[pos(23,idx,idy)]
// ical
#define         d                   state_var[pos(24,idx,idy)]
#define         ff                  state_var[pos(25,idx,idy)]
#define         fs                  state_var[pos(26,idx,idy)]
#define         fcaf                state_var[pos(27,idx,idy)]
#define         fcas                state_var[pos(28,idx,idy)]
#define         jca                 state_var[pos(29,idx,idy)]
#define         nca                 state_var[pos(30,idx,idy)]
#define         nca_i               state_var[pos(31,idx,idy)]
#define         ffp                 state_var[pos(32,idx,idy)]
#define         fcafp               state_var[pos(33,idx,idy)]
// end ical
#define         xs1                 state_var[pos(34,idx,idy)]
#define         xs2                 state_var[pos(35,idx,idy)]
#define         Jrel_np             state_var[pos(36,idx,idy)]
#define         CaMKt               state_var[pos(37,idx,idy)]
// new MM ICaL states
#define         ikr_c0              state_var[pos(38,idx,idy)]
#define         ikr_c1              state_var[pos(39,idx,idy)]
#define         ikr_c2              state_var[pos(40,idx,idy)]
#define         ikr_o               state_var[pos(41,idx,idy)]
#define         ikr_i               state_var[pos(42,idx,idy)]
#define         Jrel_p              state_var[pos(43,idx,idy)]
#define         ml     		        state_var[pos(44,idx,idy)] 		//New
#define         hl     		        state_var[pos(45,idx,idy)] 		//New


#define         dv                  state_deri[pos(1,idx,idy)]
#define         dnai                state_deri[pos(2,idx,idy)]
#define         dnass               state_deri[pos(3,idx,idy)]
#define         dki                 state_deri[pos(4,idx,idy)]
#define         dkss                state_deri[pos(5,idx,idy)]
#define         dcai                state_deri[pos(6,idx,idy)]
#define         dcass               state_deri[pos(7,idx,idy)]
#define         dcansr              state_deri[pos(8,idx,idy)]
#define         dcajsr              state_deri[pos(9,idx,idy)]

#define         dm                  state_deri[pos(10,idx,idy)]
#define         dhp                 state_deri[pos(11,idx,idy)]
#define         dh                  state_deri[pos(12,idx,idy)]
#define         dj                  state_deri[pos(13,idx,idy)]

#define         djp                 state_deri[pos(14,idx,idy)]
#define         dmL                 state_deri[pos(15,idx,idy)]
#define         dhL                 state_deri[pos(16,idx,idy)]
#define         dhLp                state_deri[pos(17,idx,idy)]
#define         da                  state_deri[pos(18,idx,idy)]
#define         diF                 state_deri[pos(19,idx,idy)]
#define         diS                 state_deri[pos(20,idx,idy)]
#define         dap                 state_deri[pos(21,idx,idy)]
#define         diFp                state_deri[pos(22,idx,idy)]
#define         diSp                state_deri[pos(23,idx,idy)]
// ical
#define         dd                  state_deri[pos(24,idx,idy)]
#define         dff                 state_deri[pos(25,idx,idy)]
#define         dfs                 state_deri[pos(26,idx,idy)]
#define         dfcaf               state_deri[pos(27,idx,idy)]
#define         dfcas               state_deri[pos(28,idx,idy)]
#define         djca                state_deri[pos(29,idx,idy)]
#define         dnca                state_deri[pos(30,idx,idy)]
#define         dnca_i              state_deri[pos(31,idx,idy)]
#define         dffp                state_deri[pos(32,idx,idy)]
#define         dfcafp              state_deri[pos(33,idx,idy)]


#define         dxs1                state_deri[pos(34,idx,idy)]
#define         dxs2                state_deri[pos(35,idx,idy)]
#define         dJrel_np            state_deri[pos(36,idx,idy)]
#define         dCaMKt              state_deri[pos(37,idx,idy)]
// new MM ICaL states
#define         dikr_c0             state_deri[pos(38,idx,idy)]
#define         dikr_c1             state_deri[pos(39,idx,idy)]
#define         dikr_c2             state_deri[pos(40,idx,idy)]
#define         dikr_o              state_deri[pos(41,idx,idy)]
#define         dikr_i              state_deri[pos(42,idx,idy)]
#define         dJrel_p             state_deri[pos(43,idx,idy)]
#define         dml     		    state_deri[pos(44,idx,idy)] 		//New
#define         dhl     		    state_deri[pos(45,idx,idy)] 		//New


////////////////////////////////////////////////////////////////////////////////
//! Compute an Euler step of the LR1 model using the Rush&Larsen method.
////////////////////////////////////////////////////////////////////////////////

/* -------------------------------------------------------------
 * Model currents
 */
__device__ double comp_current(double *state_var,double *state_deri, int idx, int idy,
					double GCAL_C, double GKS_C, double GKR_C, double GNAL_opt);

__device__ double comp_current(double *state_var,double *state_deri, int idx, int idy,
					double GCAL_C, double GKS_C, double GKR_C, double GNAL_opt)
//----------------------------------------------------------------------- Paras -----------------------------------------------------------------------------------------------
{

        //update CaMK
//!--------------------------------------------------------
        double CaMKb=(CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass));
        double CaMKa=(CaMKb+CaMKt);
        dCaMKt=(aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt);

        // INa
//! 1 ---------------------------------------------------------
        //[INa, dm, dh, dhp, dj, djp] = getINa_Grandi(v, m, h, hp, j, jp, fINap, ENa, INa_Multiplier);
        // The Grandi implementation updated with INa phosphorylation.
        double ah,bh,aj,bj;
        // m gate
        double mss = 1 / ((1 + exp( -(56.86 + v) / 9.03 ))*(1 + exp( -(56.86 + v) / 9.03 )));
        double taum = 0.1292 * exp(-((v+45.79)/15.54)*((v+45.79)/15.54)) + 0.06487 * exp(-((v-4.823)/51.12)*((v-4.823)/51.12));
        //dm = (mss - m) / taum;
        dm = (mss-(mss-m)*exp(-dt/taum)-m)/dt;

        // h gate
        if(v >= -40.0)
        {
                ah = (0.0);
                bh = (0.77 / (0.13*(1 + exp( -(v + 10.66) / 11.1 ))));
        }
        else
        {
                ah = (0.057 * exp( -(v + 80) / 6.8 ));
                bh = ((2.7 * exp( 0.079 * v) + 3.1*1e5 * exp(0.3485 * v)));
        }
        double tauh = 1 / (ah + bh);
        double hss = 1 / ((1 + exp( (v + 71.55)/7.43 ))*(1 + exp( (v + 71.55)/7.43 )));
        //dh = (hss - h) / tauh;
        dh = (hss-(hss-h)*exp(-dt/tauh)-h)/dt;

        // j gate
        if(v >= -40)
        {
                aj = (0.0) ;
                bj = ((0.6 * exp( 0.057 * v)) / (1 + exp( -0.1 * (v + 32) )));
        }
        else
        {
                aj = (((-2.5428 * 1e4*exp(0.2444*v) - 6.948*1e-6 * exp(-0.04391*v)) * (v + 37.78)) /
                (1 + exp( 0.311 * (v + 79.23) )));
                bj = ((0.02424 * exp( -0.01052 * v )) / (1 + exp( -0.1378 * (v + 40.14) )));
        }

        double tauj = 1 / (aj + bj);
        double jss = 1 / ((1 + exp( (v + 71.55)/7.43 ))*(1 + exp( (v + 71.55)/7.43 )));
        //dj = (jss - j) / tauj;
        dj = (jss-(jss-j)*exp(-dt/tauj)-j)/dt;

        // h phosphorylated
        double hssp = 1 / ((1 + exp( (v + 71.55 + 6)/7.43 ))*(1 + exp( (v + 71.55 + 6)/7.43 )));
        //dhp = (hssp - hp) / tauh;
        dhp = (hssp-(hssp-hp)*exp(-dt/tauh)-hp)/dt;

        // j phosphorylated
        double taujp = 1.46 * tauj;
        //djp = (jss - jp) / taujp;
        djp = (jss-(jss-jp)*exp(-dt/taujp)-jp)/dt;

        double GNa = 11.7802;
        double INa=GNa*(v-ENa)*m*m*m*((1.0-fINap)*h*j+fINap*hp*jp);

        //!     New adding I_NaL in HRd
		double 	gbarnal	= 0.0065;
		double  alphaml = 0.32*(v + 47.13)/(1.0 - exp(-0.1*(v + 47.13)));
		double  betaml  = 0.08*exp(-v/11.0);
		double  mltau   = 1.0/(alphaml + betaml);
		double  mlss    = alphaml/(alphaml + betaml);
		double  hlss    = 1.0/(1.0 + exp((v + 91.0)/6.1));
		double  hltau   = 600.0;
		double  enal    = ENa;
		//---------------------------------------------------------
		double xinal   = GNAL_opt*gbarnal*ml*ml*ml*hl*(v-enal);
		dml = (mlss-(mlss-ml)*exp(-dt/mltau)-ml)/dt;
		dhl = (hlss-(hlss-hl)*exp(-dt/hltau)-hl)/dt;
//!------------------------------------------------------------------------------------------ I_NaL in HRd ------------------------------------------------------------------
        INa=INa+xinal;

        //calculate INaL
//! 2 ---------------------------------------------------------
        //function [INaL,dmL,dhL,dhLp] = getINaL_ORd2011
        //(v, mL, hL, hLp, fINaLp, ENa, celltype,
           //INaL_Multiplier)
        double mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
        double tm = 0.1292 * exp(-((v+45.79)/15.54)*((v+45.79)/15.54)) + 0.06487 * exp(-((v-4.823)/51.12)*((v-4.823)/51.12));
        double tmL=tm;
        //dmL=(mLss-mL)/tmL;
        dmL = (mLss-(mLss-mL)*exp(-dt/tmL)-mL)/dt;
        double hLss=1.0/(1.0+exp((v+87.61)/7.488));
        double thL=200.0;
        //dhL=(hLss-hL)/thL;
        dhL = (hLss-(hLss-hL)*exp(-dt/thL)-hL)/dt;
        double hLssp=1.0/(1.0+exp((v+93.81)/7.488));
        double thLp=3.0*thL;
        //dhLp=(hLssp-hLp)/thLp;
        dhLp = (hLssp-(hLssp-hLp)*exp(-dt/thLp)-hLp)/dt;
        double GNaL=0.0279;
        if (celltype==1) GNaL=GNaL*0.6;

        double INaL=GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);


        // ITo
//! 3 ---------------------------------------------------------
        //function [Ito,da,diF,diS,dap,diFp, diSp] = getITo_ORd2011(v, a, iF, iS, ap, iFp, iSp, fItop, EK, celltype, Ito_Multiplier)
        //calculate Ito
        double ass=1.0/(1.0+exp((-(v-14.34))/14.82));
        double ta=1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));
        da=(ass-a)/ta;
        double iss=1.0/(1.0+exp((v+43.94)/5.711));
        double delta_epi;
        if(celltype==1)
        {
            delta_epi=1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
        }
        else
        {
            delta_epi=1.0;
        }
        double tiF=4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
        double tiS=23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
        tiF=tiF*delta_epi;
        tiS=tiS*delta_epi;
        double AiF=1.0/(1.0+exp((v-213.6)/151.2));
        double AiS=1.0-AiF;
        diF=(iss-iF)/tiF;
        diS=(iss-iS)/tiS;
        double i=AiF*iF+AiS*iS;
        double assp=1.0/(1.0+exp((-(v-24.34))/14.82));
        dap=(assp-ap)/ta;
        double dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
        double dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
        double tiFp=dti_develop*dti_recover*tiF;
        double tiSp=dti_develop*dti_recover*tiS;
        diFp=(iss-iFp)/tiFp;
        diSp=(iss-iSp)/tiSp;
        double ip=AiF*iFp+AiS*iSp;
        double Gto=0.16;
        if(celltype!=0) Gto=Gto*2.0;

        double Ito=Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);


        // ICaL
//! 4 ---------------------------------------------------------
        // a variant updated by jakub, using a changed activation curve
        // it computes both ICaL in subspace and myoplasm (_i)
        //function [ICaL_ss,ICaNa_ss,ICaK_ss,ICaL_i,ICaNa_i,ICaK_i,dd,dff,dfs,dfcaf,dfcas,...
            //djca,dnca, dnca_i, dffp,dfcafp, PhiCaL_ss, PhiCaL_i, gammaCaoMyo, gammaCaiMyo] = getICaL_ORd2011_jt(v, d,ff,fs,fcaf,fcas,jca,nca, nca_i,ffp,fcafp,...
            //fICaLp, cai, cass, cao, nai, nass, nao, ki,kss,ko, cli, clo, celltype, ICaL_fractionSS, ICaL_PCaMultiplier)

        //calculate ICaL, ICaNa, ICaK

        double dss=1.0763*exp(-1.0070*exp(-0.0829*(v)));  //magyar
        if(v >31.4978) // activation cannot be greater than 1
        {
                dss = 1.0;
        }

        double td= 0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));

        dd=(dss-d)/td;
        double fss=1.0/(1.0+exp((v+19.58)/3.696));
        double tff=7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
        double tfs=1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
        double Aff=0.6;
        double Afs=1.0-Aff;
        dff=(fss-ff)/tff;
        dfs=(fss-fs)/tfs;
        double f=Aff*ff+Afs*fs;
        double fcass=fss;
        double tfcaf=7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
        double tfcas=100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));

        double Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));

        double Afcas=1.0-Afcaf;
        dfcaf=(fcass-fcaf)/tfcaf;
        dfcas=(fcass-fcas)/tfcas;
        double fca=Afcaf*fcaf+Afcas*fcas;

        double tjca = 75;
        double jcass = 1.0/(1.0+exp((v+18.08)/(2.7916)));
        djca=(jcass-jca)/tjca;
        double tffp=2.5*tff;
        dffp=(fss-ffp)/tffp;
        double fp=Aff*ffp+Afs*fs;
        double tfcafp=2.5*tfcaf;
        dfcafp=(fcass-fcafp)/tfcafp;
        double fcap=Afcaf*fcafp+Afcas*fcas;

        // SS nca
        double Kmn=0.002;
        double k2n=500.0;
        double km2n=jca*1;
        double anca=1.0/(k2n/km2n+pow((1.0+Kmn/cass),4.0));
        dnca=anca*k2n-nca*km2n;

        // myoplasmic nca
        double anca_i = 1.0/(k2n/km2n+pow((1.0+Kmn/cai),4.0));
        dnca_i = anca_i*k2n-nca_i*km2n;

        // SS driving force
        double Io = 0.5*(nao + ko + clo + 4*cao)/1000 ; // ionic strength outside. /1000 is for things being in micromolar
        double Ii = 0.5*(nass + kss + cli + 4*cass)/1000 ; // ionic strength outside. /1000 is for things being in micromolar
        // The ionic strength is too high for basic DebHuc. We'll use Davies
        double dielConstant = 74.0; // water at 37?
        double temp = 310.0; // body temp in kelvins.
        double constA = 1.82*1e6*pow((dielConstant*temp),(-1.5));

        double gamma_cai = exp(-constA * 4.0 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
        double gamma_cao = exp(-constA * 4.0 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
        double gamma_nai = exp(-constA * 1.0 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
        double gamma_nao = exp(-constA * 1.0 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
        double gamma_ki = exp(-constA * 1.0 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
        double gamma_kao = exp(-constA * 1.0 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));


        double PhiCaL_ss =  4.0*vffrt*(gamma_cai*cass*exp(2.0*vfrt)-gamma_cao*cao)/(exp(2.0*vfrt)-1.0);
        double PhiCaNa_ss =  1.0*vffrt*(gamma_nai*nass*exp(1.0*vfrt)-gamma_nao*nao)/(exp(1.0*vfrt)-1.0);
        double PhiCaK_ss =  1.0*vffrt*(gamma_ki*kss*exp(1.0*vfrt)-gamma_kao*ko)/(exp(1.0*vfrt)-1.0);

     
        // Myo driving force
        Io = 0.5*(nao + ko + clo + 4*cao)/1000 ; // ionic strength outside. /1000 is for things being in micromolar
        Ii = 0.5*(nai + ki + cli + 4*cai)/1000 ; // ionic strength outside. /1000 is for things being in micromolar
        // The ionic strength is too high for basic DebHuc. We'll use Davies


        gamma_cai = exp(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
        gamma_cao = exp(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
        gamma_nai = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
        gamma_nao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
        gamma_ki = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
        double amma_kao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));

        double gammaCaoMyo = gamma_cao;
        double gammaCaiMyo = gamma_cai;

        double PhiCaL_i =  4.0*vffrt*(gamma_cai*cai*exp(2.0*vfrt)-gamma_cao*cao)/(exp(2.0*vfrt)-1.0);
        double PhiCaNa_i =  1.0*vffrt*(gamma_nai*nai*exp(1.0*vfrt)-gamma_nao*nao)/(exp(1.0*vfrt)-1.0);
        double PhiCaK_i =  1.0*vffrt*(gamma_ki*ki*exp(1.0*vfrt)-gamma_kao*ko)/(exp(1.0*vfrt)-1.0);

       
        // The rest
        double PCa=8.3757e-05*GCAL_C;

        if(celltype==1)
        {
            PCa=PCa*1.2;
        }
        else if(celltype==2)
        {
                PCa=PCa*2;
        }

        double PCap=1.1*PCa;
        double PCaNa=0.00125*PCa;
        double PCaK=3.574e-4*PCa;
        double PCaNap=0.00125*PCap;
        double PCaKp=3.574e-4*PCap;

        double ICaL_ss=(1.0-fICaLp)*PCa*PhiCaL_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
        double ICaNa_ss=(1.0-fICaLp)*PCaNa*PhiCaNa_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
        double ICaK_ss=(1.0-fICaLp)*PCaK*PhiCaK_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK_ss*d*(fp*(1.0-nca)+jca*fcap*nca);

        double ICaL_i=(1.0-fICaLp)*PCa*PhiCaL_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCap*PhiCaL_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
        double ICaNa_i=(1.0-fICaLp)*PCaNa*PhiCaNa_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCaNap*PhiCaNa_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
        double ICaK_i=(1.0-fICaLp)*PCaK*PhiCaK_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCaKp*PhiCaK_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);


        // And we weight ICaL (in ss) and ICaL_i
        ICaL_i = ICaL_i * (1.0-ICaL_fractionSS);
        ICaNa_i = ICaNa_i * (1.0-ICaL_fractionSS);
        ICaK_i = ICaK_i * (1.0-ICaL_fractionSS);
        ICaL_ss = ICaL_ss * ICaL_fractionSS;
        ICaNa_ss = ICaNa_ss * ICaL_fractionSS;
        ICaK_ss = ICaK_ss * ICaL_fractionSS;

        double ICaL = 1*(ICaL_ss + ICaL_i);
        double ICaNa = ICaNa_ss + ICaNa_i;
        double ICaK = ICaK_ss + ICaK_i;
        double ICaL_tot = ICaL + ICaNa + ICaK;
		
		ICaL_tot = ICaL_tot;
		
        // I_Kr
//! 5 ---------------------------------------------------------
        //Variant based on Lu-Vandenberg
        //function [IKr, dc0, dc1, dc2, do, di ] = getIKr_ORd2011_MM(V,c0,c1, c2, o, i,...
            //ko, EK, celltype, IKr_Multiplier)

        //Extracting state vector
        // c3 = y(1);
        // c2 = y(2);
        // c1 = y(3);
        // o = y(4);
        // i = y(5);
        double b = 0.0; // no channels blocked in via the mechanism of specific MM states
		b = b;

        // transition rates
        // from c0 to c1 in l-v model,
        double alpha = 0.1161 * exp(0.2990 * vfrt);
        // from c1 to c0 in l-v/
        double beta =  0.2442 * exp(-1.604 * vfrt);

        // from c1 to c2 in l-v/
        double alpha1 = 1.25 * 0.1235 ;
        // from c2 to c1 in l-v/
        double beta1 =  0.1911;

        // from c2 to o/           c1 to o
        double alpha2 =0.0578 * exp(0.9710 * vfrt); //
        // from o to c2/
        double beta2 = 0.349e-3* exp(-1.062 * vfrt); //

        // from o to i
        double alphai = 0.2533 * exp(0.5953 * vfrt); //
        // from i to o
        double betai = 1.25* 0.0522 * exp(-0.8209 * vfrt); //

        // from c2 to i (from c1 in orig)
        double alphac2ToI = 0.52e-4 * exp(1.525 * vfrt); //
        // from i to c2
        // betaItoC2 = 0.85e-8 * exp(-1.842 * vfrt); %
        double betaItoC2 = (beta2 * betai * alphac2ToI)/(alpha2 * alphai); //
        // transitions themselves
        // for reason of backward compatibility of naming of an older version of a
        // MM IKr, c3 in code is c0 in article diagram, c2 is c1, c1 is c2.

        dikr_c0 = ikr_c1 * beta - ikr_c0 * alpha; // delta for c0
        dikr_c1 = ikr_c0 * alpha + ikr_c2*beta1 - ikr_c1*(beta+alpha1); // c1
        dikr_c2 = ikr_c1 * alpha1 + ikr_o*beta2 + ikr_i*betaItoC2 - ikr_c2 * (beta1 + alpha2 + alphac2ToI);
        // subtraction is into c2, to o, to i. % c2
        dikr_o = ikr_c2 * alpha2 + ikr_i*betai - ikr_o*(beta2+alphai);
        dikr_i = ikr_c2*alphac2ToI + ikr_o*alphai - ikr_i*(betaItoC2 + betai);

        double GKr = 0.0321 * sqrt(ko/5.0)*GKR_C; // 1st element compensates for change to ko (sqrt(5/5.4)* 0.0362)
        if(celltype==1)
        {
            GKr=GKr*1.3;
        }
        else if(celltype==2)
        {
            GKr=GKr*0.8;
        }

        double IKr = GKr * ikr_o  * (v-EK);


        //calculate IKs
//! 6 ---------------------------------------------------------
        //function [IKs,dxs1, dxs2] = getIKs_ORd2011
        //(v,xs1, xs2, cai, EKs, celltype, IKs_Multiplier)
        double xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
        double txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
        dxs1=(xs1ss-xs1)/txs1;
        double xs2ss=xs1ss;
        double txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
        dxs2=(xs2ss-xs2)/txs2;
        double KsCa=1.0+0.6/(1.0+pow((3.8e-5/cai),1.4));
        double GKs= 0.0011*GKS_C;
        if (celltype==1) GKs=GKs*1.4;

        double IKs=GKs*KsCa*xs1*xs2*(v-EKs);


        //% IK1
//! 7 ---------------------------------------------------------
        //function [IK1] = getIK1_CRLP(v,  ko , EK,
        //celltype, IK1_Multiplier)
        double aK1 = 4.094/(1+exp(0.1217*(v-EK-49.934)));
        double bK1 = (15.72*exp(0.0674*(v-EK-3.257))+exp(0.0618*(v-EK-594.31)))/(1+exp(-0.1629*(v-EK+14.207)));
        double K1ss = aK1/(aK1+bK1);

        double GK1=0.6992; //0.7266; %* sqrt(5/5.4))
        if (celltype==1)
            GK1=GK1*1.2;
        else if (celltype==2)
            GK1=GK1*1.3;
        double IK1=GK1*sqrt(ko/5)*K1ss*(v-EK);


        // I_NaCa
//! 8 ---------------------------------------------------------
        //function [ INaCa_i, INaCa_ss] = getINaCa_ORd2011(v,F,R,T, nass, nai, nao, cass, cai, cao, celltype, INaCa_Multiplier, INaCa_fractionSS)
        double zca = 2.0;
        double kna1=15.0;
        double kna2=5.0;
        double kna3=88.12;
        double kasymm=12.5;
        double wna=6.0e4;
        double wca=6.0e4;
        double wnaca=5.0e3;
        double kcaon=1.5e6;
        double kcaoff=5.0e3;
        double qna=0.5224;
        double qca=0.1670;
        double hca=exp((qca*v*F)/(R*T));
        double hna=exp((qna*v*F)/(R*T));
        double h1=1+nai/kna3*(1+hna);
        double h2=(nai*hna)/(kna3*h1);
        double h3=1.0/h1;
        double h4=1.0+nai/kna1*(1+nai/kna2);
        double h5=nai*nai/(h4*kna1*kna2);
        double h6=1.0/h4;
        double h7=1.0+nao/kna3*(1.0+1.0/hna);
        double h8=nao/(kna3*hna*h7);
        double h9=1.0/h7;
        double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
        double h11=nao*nao/(h10*kna1*kna2);
        double h12=1.0/h10;
        double k1=h12*cao*kcaon;
        double k2=kcaoff;
        double k3p=h9*wca;
        double k3pp=h8*wnaca;
        double k3=k3p+k3pp;
        double k4p=h3*wca/hca;
        double k4pp=h2*wnaca;
        double k4=k4p+k4pp;
        double k5=kcaoff;
        double k6=h6*cai*kcaon;
        double k7=h5*h2*wna;
        double k8=h8*h11*wna;
        double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
        double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
        double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
        double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
        double E1=x1/(x1+x2+x3+x4);
        double E2=x2/(x1+x2+x3+x4);
        double E3=x3/(x1+x2+x3+x4);
        double E4=x4/(x1+x2+x3+x4);
        double KmCaAct=150.0e-6;
        double allo=1.0/(1.0+(KmCaAct/cai)*(KmCaAct/cai));
        double zna=1.0;
        double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
        double JncxCa=E2*k2-E1*k1;
        double Gncx= 0.0034;
        if (celltype==1)
            Gncx=Gncx*1.1;
        else if (celltype==2)
            Gncx=Gncx*1.4;
        double INaCa_i=(1-INaCa_fractionSS)*Gncx*allo*(zna*JncxNa+zca*JncxCa);

        //calculate INaCa_ss
        h1=1+nass/kna3*(1+hna);
        h2=(nass*hna)/(kna3*h1);
        h3=1.0/h1;
        h4=1.0+nass/kna1*(1+nass/kna2);
        h5=nass*nass/(h4*kna1*kna2);
        h6=1.0/h4;
        h7=1.0+nao/kna3*(1.0+1.0/hna);
        h8=nao/(kna3*hna*h7);
        h9=1.0/h7;
        h10=kasymm+1.0+nao/kna1*(1+nao/kna2);
        h11=nao*nao/(h10*kna1*kna2);
        h12=1.0/h10;
        k1=h12*cao*kcaon;
        k2=kcaoff;
        k3p=h9*wca;
        k3pp=h8*wnaca;
        k3=k3p+k3pp;
        k4p=h3*wca/hca;
        k4pp=h2*wnaca;
        k4=k4p+k4pp;
        k5=kcaoff;
        k6=h6*cass*kcaon;
        k7=h5*h2*wna;
        k8=h8*h11*wna;
        x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
        x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
        x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
        x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
        E1=x1/(x1+x2+x3+x4);
        E2=x2/(x1+x2+x3+x4);
        E3=x3/(x1+x2+x3+x4);
        E4=x4/(x1+x2+x3+x4);
        KmCaAct=150.0e-6 ;
        allo=1.0/(1.0+(KmCaAct/cass)*(KmCaAct/cass));
        JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
        JncxCa=E2*k2-E1*k1;
        double INaCa_ss=INaCa_fractionSS*Gncx*allo*(zna*JncxNa+zca*JncxCa);


        //%calculate INaK
//! 9 ---------------------------------------------------------
        //function INaK = getINaK_ORd2011(v, F, R, T,
        //nai, nao, ki, ko, celltype, INaK_Multiplier)
        zna=1.0;
        double k1p=949.5;
        double k1m=182.4;
        double k2p=687.2;
        double k2m=39.4;
        k3p=1899.0;
        double k3m=79300.0;
        k4p=639.0;
        double k4m=40.0;
        double Knai0=9.073;
        double Knao0=27.78;
        double delta=-0.1550;
        double Knai=Knai0*exp((delta*v*F)/(3.0*R*T));
        double Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
        double Kki=0.5;
        double Kko=0.3582;
        double MgADP=0.05;
        double MgATP=9.8;
        double Kmgatp=1.698e-7;
        double H=1.0e-7;
        double eP=4.2;
        double Khp=1.698e-7;
        double Knap=224.0;
        double Kxkur=292.0;
        double P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
        double a1=(k1p*pow((nai/Knai),3.0))/(pow((1.0+nai/Knai),3.0)+(1.0+ki/Kki)*(1.0+ki/Kki)-1.0);
        double b1=k1m*MgADP;
        double a2=k2p;
        double b2=(k2m*pow((nao/Knao),3.0))/((pow((1.0+nao/Knao),3.0))+(1.0+ko/Kko)*(1.0+ko/Kko)-1.0);
        double a3=(k3p*(ko/Kko)*(ko/Kko))/(pow((1.0+nao/Knao),3.0)+(1.0+ko/Kko)*(1.0+ko/Kko)-1.0);
        double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
        double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
        double b4=(k4m*(ki/Kki)*(ki/Kki))/(pow((1.0+nai/Knai),3.0)+(1.0+ki/Kki)*(1.0+ki/Kki)-1.0);
        x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
        x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
        x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
        x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
        E1=x1/(x1+x2+x3+x4);
        E2=x2/(x1+x2+x3+x4);
        E3=x3/(x1+x2+x3+x4);
        E4=x4/(x1+x2+x3+x4);
        double zk=1.0;
        double JnakNa=3.0*(E1*a3-E2*b3);
        double JnakK=2.0*(E4*b1-E3*a1);
        double Pnak= 15.4509;
        if (celltype==1)
            Pnak=Pnak*0.9;
        else if (celltype==2)
            Pnak=Pnak*0.7;

        double INaK=Pnak*(zna*JnakNa+zk*JnakK);


        // ---------------------------------------------Minor/background currents----------------------------------------------------------------------------------------------------------------------
		
        //calculate IKb
//!---------------------------------------------------------
        double xkb=1.0/(1.0+exp(-(v-10.8968)/(23.9871)));
        double GKb=0.0189;
        if (celltype==1)
            GKb=GKb*0.6;
        double IKb=GKb*xkb*(v-EK);

        //calculate INab
        double PNab=1.9239e-09;
        double INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);

        //calculate ICab
        double PCab=5.9194e-08;
        double ICab=PCab*4.0*vffrt*(gammaCaiMyo*cai*exp(2.0*vfrt)-gammaCaoMyo*cao)/(exp(2.0*vfrt)-1.0);


//!---------------------------------------------------------

        //calculate IpCa
        double GpCa=5e-04;
        double IpCa=GpCa*cai/(0.0005+cai);

        // Chloride
        // I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current

        double ecl = (R*T/F)*log(cli/clo);            // [mV]

        double Fjunc = 1.0;
        double Fsl = 1.0-Fjunc;
        // fraction in SS and in myoplasm - as per literature, I(Ca)Cl is in junctional subspace

        Fsl = 1.0-Fjunc; // fraction in SS and in myoplasm
        double GClCa =  0.2843;   // [mS/uF]
        double GClB =  1.98e-3;        // [mS/uF] %
        double KdClCa = 0.1;    // [mM]

        double I_ClCa_junc = Fjunc*GClCa/(1.0+KdClCa/cass)*(v-ecl);
        double I_ClCa_sl = Fsl*GClCa/(1.0+KdClCa/cai)*(v-ecl);


        double I_ClCa = I_ClCa_junc+I_ClCa_sl;
        double I_Clbk = GClB*(v-ecl);


        // Calcium handling
        //calculate ryanodione receptor calcium induced calcium release from the jsr
        double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));

        // Jrel
//!---------------------------------------------------------
        //[Jrel, dJrel_np, dJrel_p] = getJrel_ORd2011(Jrel_np, Jrel_p, ICaL_ss,cass, cajsr, fJrelp, celltype, Jrel_Multiplier);
        // Jrel
        //function [Jrel, dJrelnp, dJrelp] = getJrel_ORd2011(Jrelnp, Jrelp, ICaL, cass, cajsr, fJrelp, celltype, Jrel_Multiplier)
        double jsrMidpoint = 1.7;

        double bt=4.75;
        double a_rel=0.5*bt;
        double Jrel_inf=a_rel*(-ICaL)/(1.0+pow((jsrMidpoint/cajsr),8.0));
        if (celltype==2)
            Jrel_inf=Jrel_inf*1.7;
        double tau_rel=bt/(1.0+0.0123/cajsr);

        if (tau_rel<0.001)
            tau_rel=0.001;

        dJrel_np=(Jrel_inf-Jrel_np)/tau_rel;
        double btp=1.25*bt;
        double a_relp=0.5*btp;
        double Jrel_infp=a_relp*(-ICaL)/(1.0+pow((jsrMidpoint/cajsr),8.0));
        if (celltype==2)
            Jrel_infp=Jrel_infp*1.7;
        double tau_relp=btp/(1.0+0.0123/cajsr);

        if (tau_relp<0.001)
            tau_relp=0.001;

        dJrel_p=(Jrel_infp-Jrel_p)/tau_relp;

        double  Jrel=1.5378 * ((1.0-fJrelp)*Jrel_np+fJrelp*Jrel_p);

        double fJupp=(1.0/(1.0+KmCaMK/CaMKa));


//!---------------------------------------------------------
        //[Jup, Jleak] = getJup_ORd2011(cai, cansr, fJupp, celltype, Jup_Multiplier);

        // Jup
        //function [Jup, Jleak] = getJup_ORd2011(cai, cansr, fJupp, celltype, Jup_Multiplier)
        /*calculate serca pump, ca uptake flux
		% camkFactor = 2.4;
                    % gjup = 0.00696;
                    % Jupnp=Jup_Multiplier * gjup*cai/(cai+0.001);
                    % Jupp=Jup_Multiplier * camkFactor*gjup*cai/(cai + 8.2500e-04);
                    % if celltype==1
                    %     Jupnp=Jupnp*1.3;
                    %     Jupp=Jupp*1.3;
                    % end
                    %
                    %
                    % Jleak=Jup_Multiplier * 0.00629 * cansr/15.0;
                    % Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;*/

        //calculate serca pump, ca uptake flux
        double  Jupnp=0.005425*cai/(cai+0.00092);
        double  Jupp=2.75*0.005425*cai/(cai+0.00092-0.00017);
        if (celltype==1)
        {
            Jupnp=Jupnp*1.3;
            Jupp=Jupp*1.3;

        }

        double  Jleak=0.0048825*cansr/15.0;
        double  Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;

        //calculate tranlocation flux
        double  Jtr=(cansr-cajsr)/60.0;

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        // calculate the stimulus current, Istim
        //amp=stimAmp;
        //duration=stimDur;
        //if t<=duration
            //Istim=amp;
        //else
            //Istim=0.0;
        //end

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        //update the membrane voltage

        dv=0.0; //-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+
                //I_ClCa+I_Clbk) + st;
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        //calculate diffusion fluxes
        double JdiffNa=(nass-nai)/2.0;
        double JdiffK=(kss-ki)/2.0;
        double Jdiff=(cass-cai)/0.2;

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        //calcium buffer constants
        double cmdnmax= 0.05;
        if (celltype==1)
            cmdnmax=cmdnmax*1.3;
        double kmcmdn=0.00238;
        double trpnmax=0.07;
        double kmtrpn=0.0005;
        double BSRmax=0.047;
        double KmBSR = 0.00087;
        double BSLmax=1.124;
        double KmBSL = 0.0087;
        double csqnmax=10.0;
        double kmcsqn=0.8;

        //update intracellular concentrations, using buffers for cai, cass, cajsr
        dnai=-(ICaNa_i+INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo;
        dnass=-(ICaNa_ss+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa;

        dki=-(ICaK_i+Ito+IKr+IKs+IK1+IKb-st-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo;
        dkss=-(ICaK_ss)*Acap/(F*vss)-JdiffK;

        double Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow((kmcmdn+cai),2.0)+trpnmax*kmtrpn/pow((kmtrpn+cai),2.0));
        dcai=Bcai*(-(ICaL_i + IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo);

        double Bcass=1.0/(1.0+BSRmax*KmBSR/pow((KmBSR+cass),2.0)+BSLmax*KmBSL/pow((KmBSL+cass),2.0));
        dcass=Bcass*(-(ICaL_ss-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);

        dcansr=Jup-Jtr*vjsr/vnsr;

        double Bcajsr=1.0/(1.0+csqnmax*kmcsqn/pow((kmcsqn+cajsr),2.0));
        dcajsr=Bcajsr*(Jtr-Jrel);

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
		return (INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+
             I_ClCa+I_Clbk);
}

////////////////////////////////////////////////////////////////////////////////
//! Set the standard initial conditions (which are defined in macros above).
////////////////////////////////////////////////////////////////////////////////
int get_initial_state(double *state_var,double *state_deri)
{ 
		int ii,idx,idy;
	
		// Standard initial conditions.
		for(idy=0;idy<Blocknum;idy++)
		{
				for(idx=0;idx<Length;idx++) 
				{
						for(ii=0;ii<NDIM;ii++) state_deri[pos(ii,idx,idy)]=0.0;
				}
		 
				for(idx=0;idx<Length;idx++)
				{
						state_var[pos(1,idx,idy)]=-88.8691566357934;
						state_var[pos(2,idx,idy)]=12.0996647655188;
						state_var[pos(3,idx,idy)]=12.1000028563765;
						state_var[pos(4,idx,idy)]=142.412524737626;
						state_var[pos(5,idx,idy)]=142.412481425842;
						state_var[pos(6,idx,idy)]=7.45541572746214e-05;
						state_var[pos(7,idx,idy)]=6.50418928341426e-05;
						state_var[pos(8,idx,idy)]=1.53037019085812;
						state_var[pos(9,idx,idy)]=1.52803094224238;
						state_var[pos(10,idx,idy)]=0.000787657400526199;
						state_var[pos(11,idx,idy)]=0.674096901201792;
						state_var[pos(12,idx,idy)]=0.830658198588696;
						state_var[pos(13,idx,idy)]=0.830466744399495;
						state_var[pos(14,idx,idy)]=0.830093612199637;
						state_var[pos(15,idx,idy)]=0.000159670117055769;
						state_var[pos(16,idx,idy)]=0.528261721740178;
						state_var[pos(17,idx,idy)]=0.288775833197764;
						state_var[pos(18,idx,idy)]=0.000944249645410894;
						state_var[pos(19,idx,idy)]=0.999616956857814;
						state_var[pos(20,idx,idy)]=0.593680589620082;
						state_var[pos(21,idx,idy)]=0.000481107253796778;
						state_var[pos(22,idx,idy)]=0.999616964658062;
						state_var[pos(23,idx,idy)]=0.654092074678260;
						state_var[pos(24,idx,idy)]=8.86091322819384e-29;
						state_var[pos(25,idx,idy)]=0.999999992783113;
						state_var[pos(26,idx,idy)]=0.938965241412012;
						state_var[pos(27,idx,idy)]=0.999999992783179;
						state_var[pos(28,idx,idy)]=0.999900458262832;
						state_var[pos(29,idx,idy)]=0.999977476316330;
						state_var[pos(30,idx,idy)]=0.000492094765239740;
						state_var[pos(31,idx,idy)]=0.000833711885764158;
						state_var[pos(32,idx,idy)]=0.999999992566681;
						state_var[pos(33,idx,idy)]=0.999999992766279;
						state_var[pos(34,idx,idy)]=0.247156543918935;
						state_var[pos(35,idx,idy)]=0.000175017075236424;
						state_var[pos(36,idx,idy)]=3.90843796133124e-24;
						state_var[pos(37,idx,idy)]=0.0110752904836162;
						state_var[pos(38,idx,idy)]=0.998073652444028;
						state_var[pos(39,idx,idy)]=0.000844745297078649;
						state_var[pos(40,idx,idy)]=0.000698171876592920;
						state_var[pos(41,idx,idy)]=0.000370404872169913;
						state_var[pos(42,idx,idy)]=1.30239063420973e-05;
						state_var[pos(43,idx,idy)]=-1.88428892080206e-22;
						state_var[pos(44,idx,idy)]= 0.00111859;		 		//New_ml           			  
						state_var[pos(45,idx,idy)]= 0.339310414;			//New_hl           			  
				}
		}
		
		return(0);
}