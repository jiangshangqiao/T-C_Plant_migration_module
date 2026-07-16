%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Temporal initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Note %%%%%
%%%%% ksv==? needs to change, according to the stury area.
%%%%% Every PFT needs give initial values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICEtm1 = 0*MASKn;
% ICE_Dtm1 = 0*MASKn;
load('/work/FAC/FGSE/IDYST/npeleg/alps_eco_hyd/Run1/climate data/glacier_initial4.mat');
ICEtm1 = ICE1;
ICE_Dtm1 = ICE_D1;
%%%%% Previous condition
Ta_t = zeros(num_cell,24);
PAR_t = zeros(num_cell,24);
An_H_t = zeros(num_cell,cc_max,24);
An_L_t = zeros(num_cell,cc_max,24);
Rdark_H_t = zeros(num_cell,cc_max,24);
Rdark_L_t = zeros(num_cell,cc_max,24);
Tdp_H_t = zeros(num_cell,cc_max,24);
Tdp_L_t = zeros(num_cell,cc_max,24);
Psi_x_H_t = zeros(num_cell,cc_max,24);
Psi_l_H_t = zeros(num_cell,cc_max,24);
Psi_x_L_t = zeros(num_cell,cc_max,24);
Psi_l_L_t = zeros(num_cell,cc_max,24);
%%%
Tdp_t = zeros(num_cell,ms_max,24);
O_t = zeros(num_cell,ms_max,24);
V_t = zeros(num_cell,ms_max,24);
Pr_sno_t = zeros(num_cell,24);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Condition initial 1 step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Soil moisture
vi = [2.6215 2.6215 7.8646 13.1077 13.1077 13.1077 26.2153 26.2153 26.2153 26.2153 26.2153 26.2153 26.2153 26.2153];
Vtm1 = zeros(num_cell,ms_max);
for jk=1:ms_max
    Vtm1(:,jk) = vi(jk)*MASKn;
end
Vicetm1 = zeros(num_cell,ms_max);
%%%%% Hydrological part
Tstm1 = -3.6*MASKn; %5*MASKn;
Tdptm1 = -0.41*ones(num_cell,ms_max); %5*ones(num_cell,ms_max);
Tdamptm1 = -0.41*MASKn; %Tstm1;
TdpI_Htm1= 0*ones(num_cell,cc_max); %Tdamptm1*ones(1,cc_max);
TdpI_Ltm1= 0*ones(num_cell,cc_max); %Tdamptm1*ones(1,cc_max);
Ts_undertm1 = NaN*Tstm1;
Tdpsnowtm1 = 0*ones(num_cell,5);
%%%
SWEtm1 = 13.566*MASKn; %0*MASKn;
SNDtm1 = 0.04*MASKn; %0*MASKn;
snow_albedotm1 = 0.65*ones(num_cell,4) ;
e_snotm1 = 0.97*MASKn;
t_slstm1 = 0*MASKn;
rostm1 = 199.56*MASKn;
SP_wctm1 = 0.41*MASKn ;
tau_snotm1 = 0*MASKn;
In_Htm1 = 0*ones(num_cell,cc_max);
In_Ltm1 = 0*ones(num_cell,cc_max);
In_Littertm1 = 0*MASKn;
In_SWEtm1 = 0*MASKn;
In_urbtm1 = 0*MASKn;
In_rocktm1 = 0*MASKn;
EKtm1 = 0*MASKn;
WATtm1 = 0*MASKn;
IP_wctm1 = 0*MASKn;
Cicewtm1 = 0*MASKn;
FROCKtm1 = 0*MASKn;
Ticetm1 = 0*MASKn;
Oicetm1 =zeros(num_cell,ms_max);
Tdebtm1 = 0*ones(num_cell,md_max);
Ws_undertm1 = 1*MASKn;
%%%
Ccrown_t_tm1 = zeros(num_cell,cc_max);
%%%%%%%%%%
Tstm0 = Tstm1;
Slo_head = reshape(Slo_top,num_cell,1)*ones(1,ms_max); %%%
q_runon =  0*MASKn;
Q_channel = zeros(m_cell,n_cell);
Qi_in = zeros(num_cell,ms_max);
Qi_out = zeros(num_cell,ms_max);
Qi_out_Rout = zeros(m_cell,n_cell,ms_max);
Qi_in_Rout = zeros(m_cell,n_cell,ms_max);
Q_exit = 0;
Qsub_exit = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% High vegetation
Citm1_sunH = 370*ones(num_cell,cc_max);
Citm1_shdH = 370*ones(num_cell,cc_max);
LAI_Htm1 = zeros(num_cell,cc_max);
B_Htm1 =  zeros(num_cell,cc_max,8);
PHE_S_Htm1 = zeros(num_cell,cc_max);
NPP_Htm1 = zeros(num_cell,cc_max);
NPPI_Htm1 = zeros(num_cell,cc_max);
Bfac_weekHtm1 = ones(num_cell,cc_max);
dflo_Htm1 = zeros(num_cell,cc_max);
AgeL_Htm1 = zeros(num_cell,cc_max);
AgeDL_Htm1 = zeros(num_cell,cc_max);
e_rel_Htm1 = ones(num_cell,cc_max);
e_relN_Htm1 = ones(num_cell,cc_max);
SAI_Htm1 = zeros(num_cell,cc_max);
hc_Htm1 = zeros(num_cell,cc_max);
Tden_Htm1 = zeros(num_cell,cc_max);
AgePl_Htm1 = zeros(num_cell,cc_max);
Nreserve_Htm1 = 1200*ones(num_cell,cc_max);
Preserve_Htm1 = 1000*ones(num_cell,cc_max);
Kreserve_Htm1 = 1000*ones(num_cell,cc_max);
NupI_Htm1 = zeros(num_cell,cc_max,3);
FNC_Htm1 = 1*ones(num_cell,cc_max);
Vx_Htm1 = zeros(num_cell,cc_max);
Vl_Htm1 = zeros(num_cell,cc_max);
Psi_x_Htm1 = zeros(num_cell,cc_max);
Psi_l_Htm1 = zeros(num_cell,cc_max);

NuLit_Htm1 = zeros(num_cell,cc_max,3);
NBLeaf_Htm1 = zeros(num_cell,cc_max);
PARI_Htm1 = zeros(num_cell,cc_max,3);
NBLI_Htm1  = zeros(num_cell,cc_max);

TBio_H = zeros(num_cell,cc_max);

%%% Dec. forest (ksv=='A--' needs to change, according to the study area)
LAI_Htm1(ksv=='A--',1) = 0.0;
B_Htm1(ksv=='A--',1,1) = 0;
B_Htm1(ksv=='A--',1,2) = 70.9; %160.023; %286.493;    %160.023; %123.47;
B_Htm1(ksv=='A--',1,3) = 85.8; %147.620; %181.750;    %147.620; %126.47;
B_Htm1(ksv=='A--',1,4) = 48.8; %122.585; %211.892;    %122.585; %74.68;
B_Htm1(ksv=='A--',1,5) = 1.18; %1.424;      %0.524; %1.43;
B_Htm1(ksv=='A--',1,6) = 1860.8; %5723.54; %3002.964; %5959.371;   %3002.964;   %% **********
B_Htm1(ksv=='A--',1,7) = 19.1;      %0.301; %2.40;
B_Htm1(ksv=='A--',1,8) = 0;
AgeL_Htm1(ksv=='A--',1) = 0;
SAI_Htm1(ksv=='A--',1) = 0.1;
hc_Htm1(ksv=='A--',1) = 11.7;
Vx_Htm1(ksv=='A--',1) = 10;
Vl_Htm1(ksv=='A--',1) = 10;
TBio_H(ksv=='A--',1) = 20;
Ccrown_t_tm1(ksv=='A--',1) = 1;
PHE_S_Htm1(ksv=='A--',1) = 1;

% %%%% Mixed Eve Dec
% %%% if cc_max=1, this one need to be commented, because this one should be use when cc_max>=2.
% LAI_Htm1(ksv=='--C',1) = 2.8;
% B_Htm1(ksv=='--C',1,1) = 280;
% B_Htm1(ksv=='--C',1,2) = 590;
% B_Htm1(ksv=='--C',1,3) = 395;
% B_Htm1(ksv=='--C',1,4) = 395;
% AgeL_Htm1(ksv=='--C',1) = 1700;
% SAI_Htm1(ksv=='--C',1) = 0.15;
% hc_Htm1(ksv=='--C',1) = 13;
% 
% LAI_Htm1(ksv=='--C',2) = 0.0;
% B_Htm1(ksv=='--C',2,1) = 0;
% B_Htm1(ksv=='--C',2,2) = 565;
% B_Htm1(ksv=='--C',2,3) = 220;
% B_Htm1(ksv=='--C',2,4) = 410;
% AgeL_Htm1(ksv=='--C',2) = 0;
% SAI_Htm1(ksv=='--C',2) = 0.10;
% hc_Htm1(ksv=='--C',2) = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Low vegetation
Citm1_sunL = 370*ones(num_cell,cc_max);
Citm1_shdL = 370*ones(num_cell,cc_max);
LAI_Ltm1 = zeros(num_cell,cc_max);
B_Ltm1 = zeros(num_cell,cc_max,8);
PHE_S_Ltm1 = zeros(num_cell,cc_max);
NPP_Ltm1 = zeros(num_cell,cc_max);
NPPI_Ltm1 = zeros(num_cell,cc_max);
Bfac_weekLtm1 = ones(num_cell,cc_max);
dflo_Ltm1 = zeros(num_cell,cc_max);
AgeL_Ltm1 = zeros(num_cell,cc_max);
AgeDL_Ltm1 = zeros(num_cell,cc_max);
SAI_Ltm1 = zeros(num_cell,cc_max);
hc_Ltm1 = zeros(num_cell,cc_max);
Tden_Ltm1 = zeros(num_cell,cc_max);
AgePl_Ltm1 = zeros(num_cell,cc_max);
e_rel_Ltm1 = ones(num_cell,cc_max);
e_relN_Ltm1 = ones(num_cell,cc_max);
Nreserve_Ltm1 = 1200*ones(num_cell,cc_max);
Preserve_Ltm1 = 1000*ones(num_cell,cc_max);
Kreserve_Ltm1 = 1000*ones(num_cell,cc_max);
NupI_Ltm1 = zeros(num_cell,cc_max,3);
FNC_Ltm1 = 1*ones(num_cell,cc_max);
Vx_Ltm1 = zeros(num_cell,cc_max);
Vl_Ltm1 = zeros(num_cell,cc_max);
Psi_x_Ltm1 = zeros(num_cell,cc_max);
Psi_l_Ltm1 = zeros(num_cell,cc_max);

NuLit_Ltm1 = zeros(num_cell,cc_max,3);
NBLeaf_Ltm1 = zeros(num_cell,cc_max);
PARI_Ltm1 = zeros(num_cell,cc_max,3);
NBLI_Ltm1  = zeros(num_cell,cc_max);

TBio_L = zeros(num_cell,cc_max);

%%%%%% Grass
LAI_Ltm1(ksv=='--C',3) = 0.0;
B_Ltm1(ksv=='--C',3,1) = 1.521;
B_Ltm1(ksv=='--C',3,2) = 0;
B_Ltm1(ksv=='--C',3,3) = 247.09; %247.09;
B_Ltm1(ksv=='--C',3,4) = 181.88; %181.88;
B_Ltm1(ksv=='--C',3,5) = 1.026; %1.026;
B_Ltm1(ksv=='--C',3,6) = 0;
B_Ltm1(ksv=='--C',3,7) = 120.38; %120.38;
B_Ltm1(ksv=='--C',3,8) = 0;
AgeL_Ltm1(ksv=='--C',3) = 0;
SAI_Ltm1(ksv=='--C',3) = 0.001;
hc_Ltm1(ksv=='--C',3) = 0.2;
Vx_Ltm1(ksv=='--C',3) = 0;
Vl_Ltm1(ksv=='--C',3) = 100;
TBio_L(ksv=='--C',3) = 1;
Ccrown_t_tm1(ksv=='--C',3) = 1;
PHE_S_Ltm1(ksv=='--C',3) = 1;

%%% Eve. shrub
LAI_Ltm1(ksv=='-B-',2) = 2.5;
B_Ltm1(ksv=='-B-',2,1) = 120.95;
B_Ltm1(ksv=='-B-',2,2) = 210.82;
B_Ltm1(ksv=='-B-',2,3) = 146.32;
B_Ltm1(ksv=='-B-',2,4) = 178.60;
B_Ltm1(ksv=='-B-',2,5) = 1.076;
B_Ltm1(ksv=='-B-',2,6) = 1578.5; %4875.96; %1978.5; %2456.5;
B_Ltm1(ksv=='-B-',2,7) = 16.94;
B_Ltm1(ksv=='-B-',2,8) = 0;
AgeL_Ltm1(ksv=='-B-',2) = 0;
SAI_Ltm1(ksv=='-B-',2) = 0.1;
hc_Ltm1(ksv=='-B-',2) = 0.8;
Vx_Ltm1(ksv=='-B-',2) = 10;
Vl_Ltm1(ksv=='-B-',2) = 100;
TBio_L(ksv=='-B-',2) = 1;
Ccrown_t_tm1(ksv=='-B-',2) = 1;
PHE_S_Ltm1(ksv=='-B-',2) = 1;




