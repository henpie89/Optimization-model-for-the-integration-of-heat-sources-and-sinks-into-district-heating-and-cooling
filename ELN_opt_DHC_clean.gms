* created by PhD student Henrik Pieper on 09.05.2017
* from the Technical University of Denmark

$set ts "1";
$set te "8760";

Sets

p production units / air, gw, sew, sea, peak, DC /
p_nDC(p) all units except DC / air, gw, sew, sea, peak /
pns(p) no seawater and gw for source calculations / air, sew, peak /
w(p) water based heat sources / gw, sew, sea, DC /
k time steps for linearization / 1*100 /
PL time steps part load / 1*90 /
t time steps for test period in hours /%ts% * %te%/
;

********************************************************************************
* Import parameters
Parameters
cp(p,t) specific heat from all heat sources for every hour in kJ*(kg*K)^-1
d(t) heat demand for every hour in MWh*h^-1
d_s(t) sorted heat demand for every hour starting with highest in MWh*h^-1
rho(p,t) density of heat sources for every hour in kg*m3^-1
T_sea0(t) seawater temperature for every hour in K
*T_s_real(t) supply temperature for every hour in K
T_hs(p,t) heat source temperatures for every hour in K
*T_r(t) return temperature for every hour in K
V_sew(t) average volume flow rate of sewage water per day in m3*h^-1
CO2(t) Hourly CO2 emissions for DK2 from Energinet.dk in kg per MWh
p_el(t) nordpool spot market price hourly in EUR*MWh^-1 for DK2
COP_d(p) COP for design conditions
PL_load(PL) part load
d_DC(t) district cooling load
*COP(p,t) COP used
*COP_PL(p,PL) COP used
*f_hs(p) share between heat source
*eta_L(p,t) eta_L
;

* Import parameters
Execute_load 'Inputdata.gdx' cp, d=Q_heat_HP, d_s=Q_heat_HP_s, rho, T_hs, V_sew, CO2, p_el;
Execute_load 'COP_Partload.gdx' PL_load;
**Execute_load 'hs_variation.gdx' f_hs;
Execute_load 'Q_DC_demand.gdx' d_DC;

Parameter T_hs_c(p,t) temperatures heat sources for cooling;
T_hs_c(p,t) = T_hs(p,t);
Parameter rho_c(p,t) density of heat sources for every hour in kg*m3^-1;
rho_c(p,t) = rho(p,t);
Parameter cp_c(p,t) specific heat from all heat sources for every hour in kJ*(kg*K)^-1;
cp_c(p,t) = cp(p,t);

Parameter T_hs(p,t) heat source temperatures for every hour in K ;
T_hs('peak',t) = 280;
Parameter rho(p,t) heat source temperatures for every hour in K ;
rho('peak',t) = 1000;
Parameter cp(p,t) heat source temperatures for every hour in K ;
cp('peak',t) = 1;

Parameter T_hs_c(p,t) heat source temperatures for every hour in K ;
T_hs_c('peak',t) = 280;
Parameter rho_c(p,t) heat source temperatures for every hour in K ;
rho_c('peak',t) = 1000;
Parameter cp_c(p,t) heat source temperatures for every hour in K ;
cp_c('peak',t) = 1;

Scalar d_max maximum heat demand for one hour in time period in MWh*h^-1;
d_max = Smax(t,d(t));
Scalar V_sew_max maximum heat demand for one hour in time period in MWh*h^-1;
V_sew_max = Smax(t,V_sew(t));
Scalar d_c_max maximum cooling demand for one hour in time period in MWh*h^-1;
d_c_max = Smax(t,d_DC(t));

********************************************************************************
*General input data
Scalar
conv conversion of kW to MW and km to m: 1000kW=1MW / 1000 /
ex_EUR_DKK exchange rate EUR DKK average 2018 / 7.45 /
T_K conversion from degreeC in Kelvin check here if GAMS doesnt run due to division by zero /273.15 /
;

********************************************************************************
*Technical input data and heat sources
Scalar
T_gw groundwater temperature / 10 /
T_sea_min minimum seawater temperature in K / 272.65 /
;

Parameter delta_T(p) maximum temperature difference of heat source in K;
delta_T(p)=6;
delta_T('sea')=3;
delta_T('DC')=10;

Parameter delta_T_c(p) maximum temperature difference of heat source in K;
delta_T_c(p)=6;
delta_T_c('DC')=25;

********************************************************************************
* DH network
Scalar
cp_DH mean specific heat of DH water for 52.5 deg and 10 bar / 4.18 /
rho_DH mean density of DH water for 52.5 deg and 10 bar / 987.3 /
cp_DC mean specific heat of DC water for 13 deg and 10 bar / 4.188 /
rho_DC mean density of DC water for 13 deg and 10 bar / 1000 /
u_DH_max max velocity in DH pipe in m per s / 2 /
d_DH_max max diameter of DH pipe in m / 1 /
f_aux_el electricity fan power consumption factor / 0.02 /
;

Scalar T_r_min minimum DH return temperature;
T_r_min= 40 + T_K;
Scalar T_s_max maximum DH return temperature;
T_s_max= 65 + T_K;

Parameter T_s(t) supply temperature for every hour in K;
T_s(t) = T_s_max;
*T_s(t) = max(70,min(-2*(T_a(t)-T_K)+90,90))+T_K;
Parameter T_r(t) return temperature for every hour in K;
T_r(t) = T_r_min;

Parameter pipe_l(p) length of source pipe for each location in m;
pipe_l(p)=0*conv;
pipe_l('gw')=0.02*conv;
pipe_l('sew')=2*conv;
pipe_l('sea')=2.5*conv;

Parameter pipe_l_c(p) length of source pipe for each location in m;
pipe_l_c(p)=pipe_l(p);

Scalar
l_DH extra DH loss for additional piping / 0.05 /
V_gw_max estimated max volume flow rate for each borehole in m3_h / 70 /
;

*Storage
Scalar
st_c_init initial storage capacity in MWh ;
st_c_init = 0;
Scalar
st_c_init_c initial cold storage capacity in MWh ;
st_c_init_c = 0;

Scalar
st_c_max maximum storage capacity in MWh change here later or remove! / 100000 /
St_loss storage losses at each hour in MWh / 0.0001 /
st_c_max_c maximum cold storage capacity in MWh change here later or remove! / 100000 /
St_loss_c cold storage losses at each hour in MWh / 0.00001 /
;
Scalar st_char_rate charging rate of storage in MWh per h;
st_char_rate = d_max;
Scalar st_char_rate_c charging rate of storage in MWh per h;
st_char_rate_c = d_c_max;

********************************************************************************
* part load and COP
Parameter load_min(p) min volume flow rate;
load_min(p)=0;
Parameter load_max(p) min volume flow rate;
load_max(p)=1;

* part load and COP
Parameter load_min_c(p) min volume flow rate;
load_min_c(p)=0;
Parameter load_max_c(p) min volume flow rate;
load_max_c(p)=1;

Parameter T_c_i_d(p) Design inlet temperature to HP for each heat source at every hour;
T_c_i_d('sea') = 4+T_K;
T_c_i_d('sew') = 11+T_K ;
T_c_i_d('gw') = T_gw+T_K ;
T_c_i_d('air') = -12+T_K ;
T_c_i_d('peak') = 11+T_K ;
T_c_i_d('DC') = 16+T_K ;
Scalar T_r_d Design return temperature of heat sink side;
T_r_d=T_r_min;
Scalar T_s_d Design supply temperature of heat sink side;
T_s_d=T_s_max;

Parameter T_h_i_d_c(p) Design inlet temperature to chiller for each heat source at every hour;
T_h_i_d_c('air') = 35+T_K;
T_h_i_d_c('sea') = 22+T_K;
T_h_i_d_c('sew') = 25+T_K;
T_h_i_d_c('gw') = 10+T_K;
T_h_i_d_c('peak') = 100+T_K;
T_h_i_d_c('DC') = 40+T_K;
Scalar T_r_d_c Design return temperature of DC;
T_r_d_c=16+T_K;
Scalar T_s_d_c Design supply temperature of DC;
T_s_d_c=6+T_K;
Parameter T_r_c(t) return temperature for every hour in K;
T_r_c(t) = T_r_d_c;
Parameter T_s_c(t) supply temperature for every hour in K;
T_s_c(t) = T_s_d_c;

**Scalar f_hs_max / 0.604575 /;

**Parameter Q_sink_d(p) heat supply design capacity in MWh*h^-1  ;
**Q_sink_d(p)= (f_hs(p))*d_max*f_hs_max;
Parameter Q_sink_d_max(p) max heat supply capacity of heat pumps for different heat sources;
Q_sink_d_max(p)=1000;
Q_sink_d_max('gw')=5;

**Parameter Q_sink_d_c(x) cooling supply design capacity in MWh*h^-1  ;
**Q_sink_d_c(x)= (f_DC(x))*d_c_max*f_DC_max;
Parameter Q_sink_d_max_c(p) max cooling supply capacity for each unit;
Q_sink_d_max_c(p)=1000;

Parameter Q_source_d_max_c(p) max cooling supply capacity for each unit;
Q_source_d_max_c(p)=1000;
Q_source_d_max_c('gw')=5;

Parameter COP_PL(p,PL) COP used for part load;
COP_PL('sew',PL)= 0.0755*PL_Load(PL)+3.8529;
COP_PL('gw',PL)= 0.0744*PL_Load(PL)+3.795;
COP_PL('sea',PL)= 0.09475*PL_Load(PL)+3.5988;
COP_PL('air',PL)= 0.0856*PL_Load(PL)+2.7989;
COP_PL('DC',PL)= 0*PL_Load(PL)+3.964;
COP_PL('peak',PL)= 1;
Parameter COP_source(p,t) COP used for temp variation;
COP_source('sew',t)= 0.0562*(T_hs('sew',t)-T_c_i_d('sew'));
COP_source('gw',t)= 0.0553*(T_hs('gw',t)-T_c_i_d('gw'));
COP_source('sea',t)= 0.0529*(T_hs('sea',t)-T_c_i_d('sea'));
COP_source('air',t)= 0.0408*(T_hs('air',t)-T_c_i_d('air'));
COP_source('DC',t)= 0.0148*(T_hs('DC',t)-T_c_i_d('DC'));
COP_source('peak',t)= 1;
Parameter COP_sink(p,t) COP used for temp variation;
COP_sink('sew',t)= -0.029*(T_s(t)-T_s_d);
COP_sink('gw',t)= -0.0283*(T_s(t)-T_s_d);
COP_sink('sea',t)= -0.0262*(T_s(t)-T_s_d);
COP_sink('air',t)= -0.0122*(T_s(t)-T_s_d);
COP_sink('DC',t)= -0.0274*(T_s(t)-T_s_d);
COP_sink('peak',t)= 1;
Parameter COP(p,t) COP for full load;
COP(p,t)=COP_PL(p,'1') + COP_source(p,t) + COP_sink(p,t);

Parameter COP_PL_c(p,PL) COP used for part load;
COP_PL_c('air',PL)= -0.5946*PL_Load(PL)+4.8566;
COP_PL_c('sew',PL)= -0.9948*PL_Load(PL)+6.7417;
COP_PL_c('sea',PL)= -1.4082*PL_Load(PL)+7.7542;
COP_PL_c('gw',PL)= -3.5671*PL_Load(PL)+14.066;
COP_PL_c('DC',PL)= COP_PL('DC',PL)-1;
COP_PL_c('peak',PL)= 1;

Parameter COP_source_c(p,t) COP used for temp variation;
COP_source_c('air',t)$(T_hs_c('air',t)>=17+T_K)= -0.1752*(T_hs_c('air',t)-T_h_i_d_c('air'));
COP_source_c('air',t)$(T_hs_c('air',t)<17+T_K)= -0.5819*(T_hs_c('air',t)-T_h_i_d_c('air'))-0.5819*(T_h_i_d_c('air')-T_K)+17.041-COP_PL_c('air','1');
COP_source_c('sew',t) = -0.2631*(T_hs_c('sew',t)-T_h_i_d_c('sew'));
COP_source_c('sea',t)$(T_hs_c('sea',t)>=9+T_K) = -0.3191*(T_hs_c('sea',t)-T_h_i_d_c('sea'));
COP_source_c('sea',t)$(T_hs_c('sea',t)<9+T_K) = -0.6388*(T_hs_c('sea',t)-T_h_i_d_c('sea'))-0.6388*(T_h_i_d_c('sea')-T_K)+16.024-COP_PL_c('sea','1');
COP_source_c('gw',t) = -0.5062*(T_hs_c('gw',t)-T_h_i_d_c('gw'));
COP_source_c('DC',t) = COP_source('DC',t);
COP_source_c('peak',t) = 0;

Parameter COP_sink_c(p,t) COP used for temp variation;
COP_sink_c('air',t)= 0.086*(T_s_c(t)-T_s_d_c);
COP_sink_c('sew',t)= 0.1416*(T_s_c(t)-T_s_d_c);
COP_sink_c('sea',t)= 0.1744*(T_s_c(t)-T_s_d_c);
COP_sink_c('gw',t)= 0*(T_s_c(t)-T_s_d_c);
COP_sink_c('DC',t)= COP_sink('DC',t);
COP_sink_c('peak',t)= 0;

Parameter COP_c(p,t) COP for full load;
COP_c(p,t)=COP_PL_c(p,'1') + COP_source_c(p,t) + COP_sink_c(p,t);

* baseload and peakload
Scalar
*COP_base assigned value to prioritize baseload / 12 /
COP_peak assigned value to use peakload the last / 1 /
;

********************************************************************************
*electricity taxes adjust for danish conditions
Scalar
fee_trans transmission fee in DKK*MWh^-1 / 80 /
* B_high low 5.5 (44%) medium 10.68 (32%) high 16.36 (24%) = 9.764
tax_el energy tax electricity in DKK*MWh^- / 150 /
;

Parameter fee_distr(t) distribution fee in DKK*MWh^-;
fee_distr(t) = 97.64 ;
Parameter tax_c_all(t)  all taxes and tariffs in EUR*MWh^-1;
tax_c_all(t)=(fee_trans+fee_distr(t)+tax_el)/ex_EUR_DKK;

********************************************************************************

Scalar
T_y liftime of system or production units / 25 /
T_y_DH liftime of DH pipes / 30 /
T_y_st liftime of storage / 25 /
T_y_peak liftime of peak boiler / 20 /
r discount rate or rate of return / 0.04 /
;

Scalar a_c_inv used to calculate annual payments of investment;
a_c_inv = r/((1+r)*(1-(1+r)**(-T_y)));
Scalar a_c_inv_NPV used to calculate annual payments of investment;
a_c_inv_NPV = r/((1+r)*(1-(1+r)**(-T_y-1)));
Scalar a_c_inv_DH used to calculate annual payments of investment for DH pipes;
a_c_inv_DH = r/((1+r)*(1-(1+r)**(-T_y_DH)));
Scalar a_c_inv_st used to calculate annual payments of investment for storage;
a_c_inv_st = r/((1+r)*(1-(1+r)**(-T_y_st)));
Scalar a_c_inv_peak used to calculate annual payments of investment for peak boiler;
a_c_inv_peak = r/((1+r)*(1-(1+r)**(-T_y_peak)));

Scalar c_inv_stor_f investment cost factor for storage depending on size in EUR per m3;
c_inv_stor_f= 87.236 ;
Scalar c_inv_stor_min min investment cost factor for storage in EUR;
c_inv_stor_min = 204585 ;
Scalar c_inv_stor_f_c investment cost factor for storage depending on size in EUR per m3;
c_inv_stor_f_c= 1500/ex_EUR_DKK ;

Parameter m_min(p) maintenance costs for each unit in EUR per MW of installed capacity;
*m_min('base')=55366;
m_min('peak')=1020;
m_min('sea')=2000;
m_min('sew')=2000;
m_min('gw')=2000;
m_min('air')=2000;
m_min('DC')=2000;
Parameter m_f(p) maintenance costs depending on operation for each unit in EUR per MWh;
*m_f('base')=1.98;
m_f('peak')=0.5;
m_f('sea')=1.3;
m_f('sew')=1.3;
m_f('gw')=2;
m_f('air')=1;
m_f('DC')=2;
Parameter m_min_c(p) maintenance cost factor for each unit in EUR per MW of installed capacity;
m_min_c('air')=0.02;
m_min_c('sea')=0.02;
m_min_c('sew')=0.02;
m_min_c('gw')=0.03;
m_min_c('peak')=100;
m_min_c('DC')=0;

Parameter c_inv_min(p) fixed investment costs for each unit in EUR;
*c_inv_min('base')=0;
c_inv_min('peak')=0;
* HP, construction, electricity related inv,consulting, heat source
c_inv_min('sea')=0.023987+21769+150000;
c_inv_min('sew')=0.023987+21769+150000;
c_inv_min('gw')=0.023987+21769+150000;
c_inv_min('air')=0.023987+21769+150000;
c_inv_min('DC')=(925210+21769+1038600+150000)/2;

Parameter c_inv_min2(p) fixed investment costs for each unit in EUR;
*c_inv_min('base')=0;
c_inv_min2('peak')=0;
* HP, construction, electricity related inv,consulting, heat source
c_inv_min2('sea')=10850;
c_inv_min2('sew')=10850+295690;
c_inv_min2('gw')=10850+317120;
c_inv_min2('air')=10850+5.5;
c_inv_min2('DC')=0;

Parameter c_inv_min0(p) investment costs depending on capacity for each unit in EUR per MW;
c_inv_min0(p)=c_inv_min(p)+c_inv_min2(p);

Parameter c_inv_f(p) investment costs depending on capacity for each unit in EUR per MW;
*c_inv_f('base')=1800000;
c_inv_f('peak')=60000;
* HP, construction, electricity related inv, heat source
c_inv_f('sea')=336670+84311+129080;
c_inv_f('sew')=336670+84311+129080;
c_inv_f('gw')=336670+84311+129080;
c_inv_f('air')=336670+84311+129080;
c_inv_f('DC')=(359080+84311+97752)/2;

Parameter c_inv_f2(p) investment costs depending on capacity for each unit in EUR per MW;
c_inv_f2('peak')=0;
* HP, construction, electricity related inv, heat source
c_inv_f2('sea')=0;
c_inv_f2('sew')=0.93385;
c_inv_f2('gw')=89729;
c_inv_f2('air')=127380;
c_inv_f2('DC')=0;

Parameter c_inv_f0(p) investment costs depending on capacity for each unit in EUR per MW;
c_inv_f0(p)=c_inv_f(p)+c_inv_f2(p);

Parameter c_inv_min_c(p) fixed investment costs for each unit in EUR;
c_inv_min_c('air')=75000;
c_inv_min_c('sea')=75000;
c_inv_min_c('sew')=75000;
c_inv_min_c('gw')=75000;
c_inv_min_c('peak')=1000000000000000;
c_inv_min_c('DC')=(925210+21769+1038600+150000)/2;

Parameter c_inv_min_c2(p) fixed investment costs for each unit in EUR;
c_inv_min_c2('air')=75000+10850;
c_inv_min_c2('sea')=75000+10850;
c_inv_min_c2('sew')=75000+10850;
c_inv_min_c2('gw')=75000+10850;
c_inv_min_c2('peak')=1000000000000000;
c_inv_min_c2('DC')=0;

Parameter c_inv_min_c0(p) investment costs depending on capacity for each unit in EUR per MW;
c_inv_min_c0(p)=c_inv_min_c(p)+c_inv_min_c2(p);

Parameter COP_d(p) design COP;
COP_d(p)=COP_PL(p,'1');
COP_d('peak') = COP_peak;
Parameter COP_d_c(p) design COP;
COP_d_c(p)=COP_PL_c(p,'1');
COP_d('peak') = COP_peak;

Parameter c_inv_f_c(p) investment costs depending on capacity for each unit in EUR per MW;
c_inv_f_c('air')=(2000+1000)*1000/ex_EUR_DKK+129080*(COP_d('air')-1)/COP_d('air')*COP_d_c('air')/(COP_d('air')-1) ;
c_inv_f_c('sea')=(2000+1000)*1000/ex_EUR_DKK+129080*(COP_d('sea')-1)/COP_d('sea')*COP_d_c('sea')/(COP_d('sea')-1)  ;
c_inv_f_c('sew')=(2000+1000)*1000/ex_EUR_DKK+129080*(COP_d('sew')-1)/COP_d('sew')*COP_d_c('sew')/(COP_d('sew')-1)  ;
c_inv_f_c('gw')=(5000+1000)*1000/ex_EUR_DKK+129080*(COP_d('gw')-1)/COP_d('gw')*COP_d_c('gw')/(COP_d('gw')-1)  ;
c_inv_f_c('peak')=1000000000000000 ;
c_inv_f_c('DC')=(359080+84311+97752)/2;

Parameter c_inv_min_free(p) fixed investment costs for each unit in EUR;
*c_inv_min('base')=0;
c_inv_min_free('peak')=0;
* HP, construction, electricity related inv,consulting, heat source
c_inv_min_free('sea')=0;
c_inv_min_free('sew')=0+295690;
c_inv_min_free('gw')=0+317120;
c_inv_min_free('air')=0+5.5;
c_inv_min_free('DC')=0;

*Scalar c_in_DH_net investment costs in mio DKK for DH network based on HOFOR estimate / 0 /;
Scalar c_in_DH_net investment costs in mio DKK for DH network based on HOFOR estimate / 36000000 /;
Scalar C_inv_DH_net investment into DH network;
C_inv_DH_net= c_in_DH_net/ex_EUR_DKK*a_c_inv_DH/a_c_inv;

********************************************************************************
* Design conditions calculated with EES
Parameter cp_d(p) design specific heat in kJ*(kg*K)^-1;
cp_d('sea')=4.079;
cp_d('sew')=4.195;
cp_d('gw')=4.076;
cp_d('air')=1.005;
cp_d('peak')=4;
cp_d('DC')=cp_DC;
Parameter rho_d(p) design density in kg*m3^-1;
rho_d('sea')=1015;
rho_d('sew')=999.7;
rho_d('gw')=1014;
rho_d('air')=1.35;
rho_d('peak')=1000;
rho_d('DC')=rho_DC;
Scalar T_lm_h_d Design logarithmic mean temperature of heat sink side;
T_lm_h_d = (T_s_d-T_r_d)/(log(T_s_d)-log(T_r_d));
Parameter T_c_o_d(p) Design inlet temperature to HP of for each heat source at every hour;
T_c_o_d(p) = T_c_i_d(p)-delta_T(p) ;
Parameter T_lm_c_d(p) Design logarithmic mean temperature of heat source side for each heat source;
T_lm_c_d(p) = (T_c_i_d(p)-T_c_o_d(p))/(log(T_c_i_d(p))-log(T_c_o_d(p)));

Parameter cp_d_c(p) design specific heat in summer at 35 deg and 0.6 humidity in kJ*(kg*K)^-1;
cp_d_c('air')=1.046 ;
cp_d_c('sea')=4.073 ;
cp_d_c('sew')=4.181 ;
cp_d_c('gw')=4.075 ;
cp_d_c('peak')=1 ;
cp_d_c('DC')=cp_DH ;
Parameter rho_d_c(p) design density in summer at 35 deg and 0.6 humidity in kg*m3^-1;
rho_d_c('air')=1.116 ;
rho_d_c('sea')=1012 ;
rho_d_c('sew')=997.1 ;
rho_d_c('gw')=1014 ;
rho_d_c('peak')=1 ;
rho_d_c('DC')=rho_DH ;
Scalar T_lm_c_d_c Design logarithmic mean temperature of heat sink side;
T_lm_c_d_c = (T_r_d_c-T_s_d_c)/(log(T_r_d_c)-log(T_s_d_c));
Parameter T_h_o_d_c(p) Design inlet temperature to HP of for each heat source at every hour;
T_h_o_d_c(p) = T_h_i_d_c(p)+delta_T_c(p) ;
Parameter T_lm_h_d_c(p) Design logarithmic mean temperature of heat source side for each heat source;
T_lm_h_d_c(p) = (T_h_o_d_c(p)-T_h_i_d_c(p))/(log(T_h_o_d_c(p))-log(T_h_i_d_c(p)));

Parameter COP_L_d(p) Lorenz COP;
COP_L_d(p) = (T_lm_h_d/(T_lm_h_d-T_lm_c_d(p)));
COP_L_d('peak') = COP_peak;

Parameter eta_L_d(p) Lorenz efficiency;
eta_L_d(p)= COP_d(p)/COP_L_d(p);
eta_L_d('peak') = 0.2;

Parameter COP_L_d_c(p) Lorenz COP;
COP_L_d_c(p) = (T_lm_c_d_c/(T_lm_h_d_c(p)-T_lm_c_d_c));
COP_L_d('peak') = COP_peak;

Parameter eta_L_d_c(p) Lorenz efficiency;
eta_L_d_c(p)= COP_d_c(p)/COP_L_d_c(p);
eta_L_d('peak') = 0.2;

********************************************************************************
*Calculation of logarithmic mean temperatures
Parameter T_c_i(p,t) inlet temperature to HP of for each heat source at every hour;
T_c_i(p,t)=T_hs(p,t);
Parameter T_c_o(p,t) outlet temperature of HP for each heat source at every hour;
T_c_o(p,t) = T_c_i(p,t)-delta_T(p) ;
T_c_o('sea',t) = max((T_c_i('sea',t)-delta_T('sea')),T_sea_min) ;
T_c_o('sew',t) = max((T_c_i('sew',t)-delta_T('sew')),T_sea_min) ;
Parameter T_lm_h(t) logarithmic mean temperature of heat sink side;
T_lm_h(t) = (T_s(t)-T_r(t))/(log(T_s(t))-log(T_r(t)));
Parameter T_lm_c(p,t) logarithmic mean temperature of heat source side for each heat source;
T_lm_c(p,t) = (T_c_i(p,t)-T_c_o(p,t))/(log(T_c_i(p,t))-log(T_c_o(p,t)));

*Calculation of logarithmic mean temperatures
Parameter T_h_i_c(p,t) inlet temperature to HP of for each heat source at every hour;
T_h_i_c(p,t)=T_hs_c(p,t);
Parameter T_h_o_c(p,t) outlet temperature of HP for each heat source at every hour;
T_h_o_c(p,t) = T_h_i_c(p,t)+delta_T_c(p) ;
Parameter T_lm_c_c(t) logarithmic mean temperature of heat sink side;
T_lm_c_c(t) = (T_r_c(t)-T_s_c(t))/(log(T_r_c(t))-log(T_s_c(t)));
Parameter T_lm_h_c(p,t) logarithmic mean temperature of heat source side for each heat source;
T_lm_h_c(p,t) = (T_h_o_c(p,t)-T_h_i_c(p,t))/(log(T_h_o_c(p,t))-log(T_h_i_c(p,t)));

********************************************************************************

Parameter COP_L(p,t) Lorenz COP;
COP_L(p,t) = (T_lm_h(t)/(T_lm_h(t)-T_lm_c(p,t)));
COP_L('peak',t) = COP_peak;
Parameter COP(p,t) design COP;
COP('peak',t) = COP_peak;
Parameter eta_L(p,t) Lorenz efficiency;
eta_L(p,t)= COP(p,t)/COP_L(p,t);
eta_L('peak',t)=0.2;

Parameter COP_L_c(p,t) Lorenz COP;
COP_L_c(p,t) = (T_lm_c_c(t)/(T_lm_h_c(p,t)-T_lm_c_c(t)));
COP_L_c('peak',t) = COP_peak;
Parameter eta_L_c(p,t) Lorenz efficiency;
eta_L_c(p,t)= COP_c(p,t)/COP_L_c(p,t);
eta_L_c('peak',t)=0.2;

********************************************************************************
* comparison of heating and cooling demands
Parameter Q_source_DH(p,t) source capacity required to satisfy heat demand by each heat source;
Q_source_DH(p,t)=d(t)/COP(p,t)*(COP(p,t)-1);
Parameter Q_sink_DC(p,t) source capacity required to satisfy heat demand by each heat source;
Q_sink_DC(p,t)=d_DC(t)/COP_c(p,t)*(COP_c(p,t)+1);
Parameter ratio_DHC(p,t) ratio between heat source capacity and DC demand;
ratio_DHC(p,t)=Q_source_DH(p,t)/Q_sink_DC(p,t);

Parameter T_hs_c2(p,t) new heat source temperature after the use of DH;
T_hs_c2(p,t)= T_c_o(p,t);

Parameter COP_source_c2(p,t) COP used for temp variation;
COP_source_c2('air',t)$(T_hs_c('air',t)>=18+T_K)= -0.1752*(T_hs_c2('air',t)-T_h_i_d_c('air'));
COP_source_c2('air',t)$(T_hs_c('air',t)<18+T_K)= -0.5819*(T_hs_c2('air',t)-T_h_i_d_c('air'))+10.104-17.041;
COP_source_c2('sew',t) = -0.2631*(T_hs_c2('sew',t)-T_h_i_d_c('sew'));
COP_source_c2('sea',t)$(T_hs_c('sea',t)>=9+T_K) = -0.3191*(T_hs_c2('sea',t)-T_h_i_d_c('sea'));
COP_source_c2('sea',t)$(T_hs_c('sea',t)<9+T_K) = -0.6388*(T_hs_c2('sea',t)-T_h_i_d_c('sea'))+13.137-16.024;
COP_source_c2('gw',t) = -0.5062*(T_hs_c2('gw',t)-T_h_i_d_c('gw'));
COP_source_c2('peak',t) = 0;
COP_source_c2('DC',t) = 0;
Parameter COP_c2(p,t) COP for full load;
COP_c2(p,t)=COP_PL_c(p,'1') + COP_source_c2(p,t) + COP_sink_c(p,t);


********************************************************************************
* Design condition equations

**Parameter Q_source_d(p) design heat capacity of source in MWh*h^-1;
**Q_source_d(p) = Q_sink_d(p)/COP_d(p)*(COP_d(p)-1);
**Parameter m_d_source(p) design mass flow rate in kg*s^-1  ;
**m_d_source(p) = Q_source_d(p)/cp_d(p)/(T_c_i_d(p)-T_c_o_d(p))*conv;
**Parameter V_d_source(p) design volume flow rate in m3 per h;
**V_d_source(p) = m_d_source(p)/rho_d(p)*3600;

$ontext
Parameter n_bh(p) number of boreholes or pipes;
n_bh(p)=1;
n_bh('gw')=round(V_d_source('gw')/V_gw_max);
Parameter V_bh(p) max volume flow rate per borehole or pipe;
V_bh(p)$(n_bh(p)>0)=V_d_source(p)/n_bh(p);
V_bh(p)$(n_bh(p)=0)=0;

Parameter Q_source_d_c(x) design heat capacity of source in MWh*h^-1;
* switched -1 to +1 for cooling
Q_source_d_c(x) = Q_sink_d_c(x)/COP_d_c(x)*(COP_d_c(x)+1);
Parameter m_d_source_c(x) design mass flow rate in kg*s^-1  ;
* temperatures switched inlet outlet and source and sink
m_d_source_c(x) = Q_source_d_c(x)/cp_d_c(x)/(T_h_o_d_c(x)-T_h_i_d_c(x))*conv;
Parameter V_d_source_c(x) design volume flow rate in m3 per h;
V_d_source_c(x) = m_d_source_c(x)/rho_d_c(x)*3600;

Parameter n_bh_c(x) number of boreholes or pipes;
n_bh_c(x)=1;
*n_bh_c('gw')=round(V_d_source_c('gw')/V_gw_max_c);
Parameter V_bh_c(x) max volume flow rate per borehole or pipe;
V_bh_c(x)$(n_bh_c(x)>0)=V_d_source_c(x)/n_bh_c(x);
V_bh_c(x)$(n_bh_c(x)=0)=0;


*V_bh('gw')=46;
Parameter d_source_d(w) diameter heat source flow;
d_source_d(w) = sqrt(V_bh(w)/3600/u_DH_max/pi*4);
Parameter d_source_d_sew diameter heat source flow;
d_source_d_sew = sqrt(V_sew_max/3600/u_DH_max/pi*4);
$offtext

Parameter V_d_source_max(p,t) max volume flow rate for heat source in m3 per hour;
V_d_source_max(p,t)=1000000000;
V_d_source_max('sew',t)=V_sew(t);

*volume_flow_rate_sink(p,t) .. V_sink(p,t) =e= Q_hs(p,t)/(rho_DH/3600*cp_DH*(T_s(t)-T_r(t))/conv);
*Scalar d_DH_max max diameter of DH pipe in m / 1 /  ;
Parameter V_DH_max max volume flow rate in DH pipes in m3 per h;
*V_DH_max(p) = Q_sink_d(p)/(rho_DH/3600*cp_DH*(T_s_d-T_r_d)/conv);
V_DH_max = u_DH_max*pi/4*(d_DH_max)**2*3600;
**Parameter d_DH_max(p) DH pipe diameter;
**d_DH_max(p) = sqrt(V_DH_max(p)/3600/u_DH_max/pi*4);

Parameter V_d_source_max_c(p,t) max volume flow rate for heat source in m3 per hour;
V_d_source_max_c(p,t)=1000000000;
Parameter V_DH_max_c max volume flow rate in DH pipes in m3 per h;
*V_DH_max_c(x) = Q_sink_d_c(x)/(rho_DC/3600*cp_DC*(T_r_d_c-T_s_d_c)/conv);
V_DH_max_c = u_DH_max*pi/4*(d_DH_max)**2*3600;
*Parameter d_DH_max_c(x) DH pipe diameter;
*d_DH_max_c(x) = sqrt(V_DH_max_c(x)/3600/u_DH_max/pi*4);

*$offtext

********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************

$ontext
* Aux electricity consumption
Scalar
g gravity / 9.81 /
f_aux_el electricity consumption factor / 0.02 /
f_pl power factor for part load should be 2point8 / 3 /
conv_kPa_m conversion from kPa to m / 0.102/
conv_kPa_bar conversion from kPa to bar / 0.01/
conv_bar_m conversion from bar to m / 10.2/
p_HP pressure in bar of heat source in evaporator of HP / 1.5 /
tr_p_gw transmissivity efficiency when pumping NOT INCLUDED YET / 0.85 /
nu_w kinematic viscosity of water / 0.00000131 /
K_entr pipe entrance bellmouth / 0.05 /
K_90 bend 90deg / 0.75 /
K_45 bend 45deg / 0.3 /
K_v butterfly valve fully open / 0.3 /
K_no no return valve / 1 /
K_out bellmouth outlet / 0.2 /
*V_gw_1 Volume flow rate gw well / 46 /
*V_gw_2 Volume flow rate gw well / 56.3 /
*d_gw_fh diameter well flex heat / 0.125 /
;

*Parameter u_gw_1 velocity gw well;
*u_gw_1= V_gw_1/(pi/4*d_gw_fh**2)/3600;
*Parameter u_gw_2 velocity gw well;
*u_gw_2= V_gw_2/(pi/4*d_gw_fh**2)/3600;

Parameter H_st(w) static head_height difference between pump and HP;
H_st('gw')= 20;
H_st('sew')= 0;
H_st('sea')= 0;

*h_dyn_gw_max dynamic lowering of groundwater level due to pumping volume flow rate / 10 /
Parameter p_source(w) pressure at source origin;
p_source(w)=1;
p_source('gw')=1;
Parameter eta_pump(w) pump efficiency;
eta_pump(w)=0.70;
Parameter A_cross_source(w) cross sectional area of source pipe;
A_cross_source(w)=pi/4*d_source_d(w)**2;
A_cross_source('sew')=pi/4*d_DH_max('sew')**2;
Parameter L_pipe_h(p) horizontal length of pipe for heat source to HP in m;
L_pipe_h(w)=pipe_l(w);
Parameter L_pipe_tot(w) total pipe length for heat source to HP in m;
L_pipe_tot(w)= L_pipe_h(w) + H_st(w);
Parameter k_r(p) roughness of pipe (steel and PVC) in m;
k_r(w)= 0.00015;
Parameter Reynolds(w) reynolds number in pipe;
Reynolds(w)$(Q_sink_d(w)>0)=u_DH_max*d_source_d(w)/nu_w;
Reynolds('sew')$(Q_sink_d('sew')>0)=u_DH_max*d_DH_max('sew')/nu_w;
Reynolds(w)$(Q_sink_d(w)=0)=0;
Parameter test(w) test;
test(w)$(Q_sink_d(w)>0) = log(k_r(w)/3.7/d_source_d(w)+5.74/(Reynolds(w)**0.9));
test('sew')$(Q_sink_d('sew')>0) = log(k_r('sew')/3.7/d_DH_max('sew')+5.74/(Reynolds('sew')**0.9));
test(w)$(Q_sink_d(w)=0) = 0;
Parameter friction(p) friction factor in pipe;
friction(w)$(Q_sink_d(w)>0)= 0.25/power(test(w),2);
friction(w)$(Q_sink_d(w)=0) = 0;
Parameter K_p(w) loss coefficient due to pipe;
K_p(w)$(Q_sink_d(w)>0)= friction(w)*L_pipe_tot(w)/d_source_d(w);
K_p('sew')$(Q_sink_d('sew')>0)= friction('sew')*L_pipe_tot('sew')/d_DH_max('sew');
K_p(w)$(Q_sink_d(w)=0)= 0;
Parameter K_tot(w) total loss coefficient ;
K_tot('gw')= K_p('gw') + K_no +4*K_90 + K_v;
K_tot('sew')= K_p('sew') + K_no +8*K_90 + 4*K_45 + K_v;
K_tot('sea')= K_p('sea') + K_no +6*K_90 + K_v;
Parameter H_dyn(w) dynamic head;
H_dyn(w)=K_tot(w)*u_DH_max**2/2/g;
Parameter H_tot(w) total head;
H_tot(w)=H_st(w)+H_dyn(w)+(p_HP-p_source(w))*conv_bar_m;

Parameter H_st_r(w) static head_height difference between pump and HP;
H_st_r('gw')= 0;
H_st_r('sew')= 0;
H_st_r('sea')= 0;
Parameter L_pipe_h_r(p) horizontal length of pipe for heat source to HP in m;
L_pipe_h_r('sew')=L_pipe_h('sew');
L_pipe_h_r('gw')=500;
L_pipe_h_r('sea')=500;
Parameter L_pipe_tot_r(w) total pipe length for heat source to HP in m;
L_pipe_tot_r(w)= L_pipe_h_r(w) + H_st_r(w);
Parameter K_p_r(w) loss coefficient due to pipe;
K_p_r(w)$(Q_sink_d(w)>0)= friction(w)*L_pipe_tot_r(w)/d_source_d(w);
K_p_r('sew')$(Q_sink_d('sew')>0)= friction('sew')*L_pipe_tot_r('sew')/d_DH_max('sew');
K_p_r(w)$(Q_sink_d(w)=0)= 0;
Parameter K_tot_r(w) total loss coefficient ;
K_tot_r('gw')= K_p_r('gw') + K_no +2*K_90 ;
K_tot_r('sew')= K_p_r('sew') + K_no +8*K_90 + 4*K_45+ K_v ;
K_tot_r('sea')= K_p_r('sea') + K_no + 2*K_90 ;
Parameter H_dyn_r(w) dynamic head;
H_dyn_r(w)=K_tot_r(w)*u_DH_max**2/2/g;
Parameter H_tot_r(w) total head;
H_tot_r(w)=H_st_r(w)+H_dyn_r(w);

Parameter P_aux_el_r(p) aux electricity consumption after evaporator of HP;
P_aux_el_r(w) = n_bh(w)*V_bh(w)*H_tot_r(w)*g*rho_d(w)/eta_pump(w)/3600/conv/conv;
P_aux_el_r('sew') = V_DH_max('sew')*H_tot_r('sew')*g*rho_DH/eta_pump('sew')/3600/conv/conv;
P_aux_el_r('peak')=0;
P_aux_el_r('air')=0;

Parameter P_aux_el_s(p) aux electricity consumption before evaporator of HP;
*P_aux_el_d(p)=f_aux_el*Q_source_d(p);
P_aux_el_s('air')=f_aux_el*Q_source_d('air');
P_aux_el_s('peak')=0;
P_aux_el_s(w)= n_bh(w)*V_bh(w)*H_tot(w)*g*rho_d(w)/eta_pump(w)/3600/conv/conv;
P_aux_el_s('sew')= V_DH_max('sew')*H_tot('sew')*g*rho_d('sew')/eta_pump('sew')/3600/conv/conv;

Parameter P_aux_el_s_gw aux electricity consumption before evaporator of HP;
P_aux_el_s_gw$(n_bh('gw')>0)=P_aux_el_s('gw')*conv/n_bh('gw');
P_aux_el_s_gw$(n_bh('gw')=0)=0;
Parameter P_aux_el_d(p) aux electricity consumption at design conditions for full load;
P_aux_el_d(p) = P_aux_el_s(p) + P_aux_el_r(p);



* requires update when also other sources are considered ('w')
Parameter P_aux_el_s_c(x) aux electricity consumption before evaporator of HP;
P_aux_el_s_c(x)=f_aux_el*Q_source_d_c(x);

Parameter P_aux_el_d_c(x) aux electricity consumption at design conditions for full load;
*P_aux_el_d(p) = P_aux_el_s(p) + P_aux_el_r(p);
P_aux_el_d_c(x) = P_aux_el_s_c(x);

$offtext
********************************************************************************
********************************************************************************
********************************************************************************

Variables
u(p) binary variable if HP capacity above zero
u_st binary variable for storage
u_c(p) binary variable if HP capacity above zero
u_free(p) g
u_st_c binary variable for storage
*u_min(p,t) min HP load
u_min_c(p,t) min DC production unit load
*u_free_pipe_x(p) not accounting for piping twice
*u_free_pipe(p) not accounting for piping twice

Q_hs(p,t) produced heat for each production unit and hour in MWh*h^-1
P_hs(p,t) power consumption for each production unit and hour in MWh*h^-1
Q_tot(p) Total heat supply in GWh
P_tot(p) Total power consumption in GWh
St_c storage capacity
st_V storage volume in m3
St_lev(t) storage level at each hour in MWh
St_char(t) charging of storage at each hour in MWh
St_dis(t) discharging of storage at each hour in MWh
Z Total production costs in EUR
z_p(p) production costs for each HP
P_max(p) d

Q_sink_d_check(p) check HP capacity
Q_hs_check(p) check HP capacity
Q_sink_d_check_c(p) check HP capacity
Q_hs_check_c(p) check HP capacity

Q_sink_d(p) heat supply design capacity in MWh*h^-1
Q_source_d(p) design heat capacity of source in MWh*h^-1
m_d_source(p) design mass flow rate in kg*s^-1
V_d_source(p) design volume flow rate in m3 per h
Q_sink_d_c(p) heat supply design capacity in MWh*h^-1
Q_source_d_c(p) design heat capacity of source in MWh*h^-1
m_d_source_c(p) design mass flow rate in kg*s^-1
V_d_source_c(p) design volume flow rate in m3 per h

V_d_source_V_max max volume flow rate for air at Vao

Q_source(p,t) required heat capacity in MWh*h^-1
Q_sink(p,t) heat supply capacity of each production unit in MWh*h^-1
Q_source_real(p,t) real heat source
V_source(p,t) real volume flow rate
V_sink(p,t) volume flow rate DH supply line
V_sink_max(p) max volume flow rate
*d_DH(p) heat source pipe diameter

********************************************************************************
*economics
C_inv_tot(p) investment costs for each production unit in EUR
C_m(p) maintenance costs for each production unit in EUR
C_inv_stor investment costs for storage in EUR
C_DH(p) investment costs for installation of DH pipes in EUR
C_DH_pipe total investment costs for DH pipes
C_el_savings savings in distribution costs for electricity

C_el_savings_p(p) savings in distribution costs for electricity for each HP
C_el_p(p) el prod costs for each HP

C_el_tot Total_costs_el_prod
C_inv_tot_all Total_costs_investments
C_m_tot Total_costs_maintenance

* aux el
P_aux_el_t(p,t) aux el power consumption for each HP and hour
Power_aux_tot_p(p) annual paux power consumption
C_el_p_aux0(p,t) C_el_p_aux0(p)
C_el_p_aux(p) aux el costs for each heat source
C_el_aux aux el tot
load_hs(p,t) load heat source
load_hs_max(p) max load heat source
Pk(p,t,k) for linearization

* Cooling
Q_DC(p,t) DC production
St_char_c(t) cold storage charging
St_dis_c(t) cold storage discharging
St_c_c cold storage capacity
St_V_c cold storage volume
St_lev_c(t) cold storage level

Q_sink_c(p,t)  DC production unit capacity based on operation
Q_source_c(p,t) required cooling capacity in MWh*h^-1
Q_source_real_c(p,t) real heat source
V_source_c(p,t) real volume flow rate
V_sink_c(p,t) volume flow rate DH supply line
V_sink_max_c(p) max volume flow rate
P_hs_c(p,t) power consumption for each production unit and hour in MWh*h^-1
P_tot_c(p) Total power consumption in GWh
Q_tot_c(p) Total heat supply in GWh
P_max_c(p) d

V_source_c(p,t) real volume flow rate
P_aux_el_t_c(p,t) aux el power consumption for each HP and hour
Power_aux_tot_p_c(p) annual paux power consumption
load_hs_c(p,t) load heat source
load_hs_max_c(p) max load heat source
Pk_c(p,t,k) for linearization

C_inv_stor_c investment costs for storage in EUR
C_inv_tot_c(p) investment costs for each production unit in EUR
C_DH_c(p) investment costs for installation of DH pipes in EUR
C_m_c(p) maintenance costs for each production unit in EUR
C_el_p_c(p) el prod costs for each HP
C_el_tot_c  Total_costs_el_prod
C_inv_tot_all_c Total_costs_investments
C_m_tot_c  Total_costs_maintenance
C_DH_pipe_c total investment costs for DH pipes
C_el_p_aux_c(p) aux el costs for each heat source
C_el_aux_c aux el costs for each heat source

Z_p_c(p) total costs for each cooling unit

Q_free_d(p) design capacity free cooler
Q_free(p,t) hourly cooling supply free cooler
C_inv_free(p) investment costs for free cooler
C_inv_free_tot investment costs for free cooler
V_free(p,t) hourly volume flow rate free cooler
V_free_d(p) design colume flow rate
m_free(p,t) hourly mass flow rate
m_free_d(p) design mass flow rate
P_free(p,t) hourly power consumption free cooler fans
P_free_d(p) design power consumption free cooler fans
P_free_tot(p) total power cons free cooler
Q_free_tot(p) total cooling supply
*Q_free_x(p,t) hourly cooling supply free cooler and chiller

u_z_free(p)
*St_char_c_DC(t)
*St_dis_c_DC(t)
*St_char_DC(t)
*St_dis_DC(t)
;

************
Binary variable u, u_st, u_c, u_st_c;
Binary variable u_min, u_min_c,  u_free, u_free_pipe, u_free_pipe_x;
Positive Variable Q_hs, Q_sink_d_check,Q_hs_check, P_hs, Q_tot, P_tot, St_V, St_c, St_lev, St_char, St_dis, Q_sink, Q_source,
Q_source_real, V_source, V_sink, V_sink_max;
Positive Variable C_inv_tot, C_m, C_inv_stor ,C_el_tot,C_inv_tot_all,C_m_tot , C_DH, C_DH_pipe, C_el_savings_p, C_el_p;
Positive Variable z_p, P_max,P_sum_real, P_max_real_1,P_max_real_2,P_max_real_3,P_max_real_4,P_max_real_5,P_max_real_6,P_max_real_7,P_max_real_8,P_max_real_9,P_max_real_10,P_max_real_11,P_max_real_12;
Positive Variable Pk,P_aux_el_t,C_el_p_aux0,C_el_p_aux,C_el_aux, Power_aux_tot_p, load_hs_max, load_hs;
Positive Variable Z_p_c, Q_DC,St_lev_c,St_char_c,St_dis_c,St_c_c,St_V_c,Q_sink_c,Q_source_c, Q_sink_d_check_c,Q_hs_check_c, P_hs_c, P_tot_c, Q_tot_c, P_max_c, Pk_c, load_hs_c,load_hs_max_c,P_aux_el_t_c,
V_source_c,Q_source_real_c,V_source_c,V_sink_c,V_sink_max_c,C_inv_stor_c,C_inv_tot_c,C_DH_c,C_m_c,C_el_p_c,C_el_tot_c,C_inv_tot_all_c,C_m_tot_c,C_DH_pipe_c, C_el_p_aux_c,C_el_aux_c   ;
Positive Variable Q_sink_d,Q_source_d,m_d_source,V_d_source,Q_sink_d_c,Q_source_d_c,m_d_source_c,V_d_source_c;
Positive Variable Q_free_d,Q_free,C_inv_free,C_inv_free_tot,V_free,V_free_d,m_free,m_free_d,P_free,P_free_d,P_free_tot,Q_free_tot,u_z_free;
********************************************************************************
********************************************************************************
********************************************************************************

Equations
cost define objective function
cost_p(p) prod cost of each HP

********************************************************************************

* production
power_eq1(p,t)
power_eq(p,t) power consumption for each production unit and hour in MWh*h^-1
supply(p,t) heat supply of each production unit
supply_max(p,t) heat supply limit of each production unit
supply_max1(p,t)  heat supply limit of each production unit
demand(t) heat demand for every hour
supply_min(p,t) min HP load

power_tot(p) total power consumption in GWh
supply_tot(p) total heat supply in GWh
power_max(p,t) maximum electrical power

power_eq_c1(p,t)
power_eq_c(p,t) power consumption for each production unit and hour in MWh*h^-1
supply_c(p,t)
supply_max_c(p,t)
supply_max1_c(p,t)
*supply_min_c(p,t)
demand_c(t) cooling demand for every hour
power_tot_c(p) total power consumption in GWh
supply_tot_c(p) total heat supply in GWh
power_max_c(p,t) maximum electrical pow

********************************************************************************
*Storage equations
*charge(t) charging of storage
*discharge(t) discharging of storage
Max_storage(t) maximum storage capacity
storage(t) storage level for every hour
storage_init(t) storage level for 1. hour
storage_end(t) storage level for last hour
storage_cap(t) storage capacity
storage_vol storage volume

Max_storage_c(t) maximum storage capacity for cold storage
storage_cap_c(t) storage capacity
storage_init_c(t) storage level for 1. hour
storage_c(t) storage level for every hour
storage_end_c(t) storage level for last hour
storage_vol_c storage volume

********************************************************************************
* Design conditions
d_supply_cap(p) heat supply design capacity in MWh*h^-1
d_supply_cap_max(p) max heat supply design capacity in MWh*h^-1
d_source_cap(p) design heat capacity of source in MWh*h^-1
d_mass_flow_rate_source(p) design mass flow rate in kg*s^-1
d_volume_flow_rate_source(p,t) design volume flow rate in m3 per h

d_supply_cap_c_peak(p) no peak load for cooling
d_supply_cap_c(p) heat supply design capacity in MWh*h^-1
d_supply_cap_max_c(p) max heat supply design capacity in MWh*h^-1
d_source_cap_c(p) design heat capacity of source in MWh*h^-1
d_source_cap_c1(p) design heat capacity of source in MWh*h^-1
d_mass_flow_rate_source_c(p) design mass flow rate in kg*s^-1
d_volume_flow_rate_source_c(p,t) design volume flow rate in m3 per h

********************************************************************************
*Calculation of heat source and heat sink capacities
source_cap(p,t) required heat capacity in MWh*h^-1
supply_cap(p,t) heat supply capacity of each production unit in MWh*h^-1
supply_cap_peak(p,t) peak heat supply capacity of each production unit in MWh*h^-1
*supply_cap_base(p,t) base load heat supply capacity of each production unit in MWh*h^-1
*supply_cap_peak(p,t) peak load heat supply capacity of each production unit in MWh*h^-1

source_cap_c(p,t) required cooling capacity in MWh*h^-1
supply_cap_c(p,t) cooling supply capacity of each production unit in MWh*h^-1

********************************************************************************
* volume flow rates
real_source_cap(p,t) for each heat source
volume_flow_rate_source(p,t) for each heat source
volume_flow_rate_source_max(p,t) maximum available volume flow rate for each heat source in m3 per h
volume_flow_rate_sink(p,t) for DH network supply
volume_flow_rate_sink_max(p,t) for DH network supply
volume_flow_rate_sink_max_real(p,t) for DH network supply
*Diameter_source(p) heat source piping diameter
volume_flow_rate_source_max_c(p,t) maximum available volume flow rate for each heat source in m3 per h
volume_flow_rate_sink_max_c(p,t) for DH network supply

********************************************************************************
*investment costs
*inv_tot_base(p) total investment costs
*inv_tot_peak(p) total investment costs
inv_tot_hs(p) total investment costs for seawater system
inv_stor total investment costs for storage
inv_DH(pns) investment Costs for installation of DH pipes in EUR for each HP
inv_DH_sea(p) investment Costs for installation of DH pipes in EUR for each HP
inv_DH_gw(p) investment Costs for installation of DH pipes in EUR for each HP

********************************************************************************
*maintenance costs
*maintenance_base(p) maintenance costs for base load unit
*maintenance_peak(p) maintenance costs for peak load unit
maintenance_hs(p) maintenance costs for seawater unit
maintenance_hs1(p)

********************************************************************************
*total costs split into parts
Total_annual_costs_el_prod Total_costs_el_prod
Total_costs_investments Total_costs_investments
Total_annual_costs_maintenance Total_costs_maintenance
Total_investment_costs_DH_pipes DH piping costs

* Cost split into heat source and parts
annual_costs_el_prod(p)  annual costs for el prod for each HP

check_Q_sink(p,t) check HP capacity
check_Q_hs(p,t) check HP capacity
check_Q_sink_c(p,t) check HP capacity
check_Q_hs_c(p,t) check HP capacity

real_source_cap_c(p,t) for each heat source
volume_flow_rate_source_c(p,t) for each heat source
volume_flow_rate_sink_c(p,t) for DH network supply
volume_flow_rate_sink_max_real_c(p,t) for DH network supply

inv_tot_hs_c(p) total investment costs for seawater system
inv_stor_c total investment costs for storage
inv_DH_c(pns) investment Costs for installation of DH pipes in EUR for each HP
inv_DH_gw_c(p) investment Costs for installation of DH pipes in EUR for each HP
inv_DH_sea_c(p) investment Costs for installation of DH pipes in EUR for each HP
maintenance_hs_c(p) maintenance costs for seawater unit
maintenance_hs_c1(p)

annual_costs_el_prod_c(p)  annual costs for el prod for each HP
Total_annual_costs_el_prod_c Total_costs_el_prod
Total_costs_investments_c Total_costs_investments
Total_annual_costs_maintenance_c Total_costs_maintenance
Total_investment_costs_DH_pipes_c DH piping costs

cost_p_c(p) cost function for each cooling unit

*inv_free(pns) investment costs for free cooling
*inv_free1(p) investment costs for free cooling
inv_free2(p) investment costs for free cooling
*inv_free3(p)
*inv_freex(pns)
*inv_free1x(p) investment costs for free cooling
*inv_free2x(p) investment costs for free cooling
inv_free_tot investment costs for free cooling
Free_power_d(p) free cooling power design cond
Free_volume_flow_rate_d(p) vol flow free cool
Free_capacity_d(p,t)  design cap free cool
Free_mass_flow_rate_d(p)  design mass flow cool
Free_mass_flow_rate(p,t)   mass flow cool
Free_volume_flow_rate(p,t) vol free cool
Free_power(p,t) free cooling power hourly
Free_power_tot(p) total power cons free cool
Free_cooling_tot(p) total cooling supply
Free_cooling_x(p,t) total cooling supply
Chiller_cooling(p,t) eq for free cooling and chiller
Chiller_cooling1(p,t) eq for chiller
Free_mass_Flow_Rate_max(p,t) mass flow rate limit
Free_cooling_x2(p,t)
Free_capacity_d2(p,t) f
*Inv_free_piping(p) make sure piping is not accounted for twice
*Free_cooling_x3(p,t)

*DC_cond_1(p,t) conditions for DHC HP
*DC_cond_2(p,t) conditions for DHC HP
DC_cond_3(p,t) conditions for DHC HP
DC_cond_4(p,t) conditions for DHC HP
*DC_cond_5(p,t) conditions for DHC HP

u_z_free_Eq1(p)
u_z_free_Eq2(p)
u_z_free_Eq3(p)

;

********************************************************************************
********************************************************************************
********************************************************************************

* Design condition equations
d_supply_cap(p) .. Q_sink_d(p) =l= d_max;
* change here again for optimization from =e= to =l=
d_supply_cap_max(p) .. Q_sink_d(p) =l= Q_sink_d_max(p)*u(p);
d_source_cap(p) .. Q_source_d(p) =e= Q_sink_d(p)/COP_d(p)*(COP_d(p)-1);
d_mass_flow_rate_source(p) .. m_d_source(p) =e= Q_source_d(p)/cp_d(p)/(T_c_i_d(p)-T_c_o_d(p))*conv;
d_volume_flow_rate_source(p,t) .. V_d_source(p) =e= m_d_source(p)/rho_d(p)*3600;

d_supply_cap_c_peak(p) .. Q_sink_d_c('peak') =e= 0;
d_supply_cap_c(p) .. Q_sink_d_c(p) =l= d_c_max;
* change here again for optimization from =e= to =l=
d_supply_cap_max_c(p) .. Q_sink_d_c(p) =l= Q_sink_d_max_c(p)*u_c(p);
d_source_cap_c1(p) .. Q_source_d_c(p) =l= Q_source_d_max_c(p)*u_c(p);
d_source_cap_c(p) .. Q_source_d_c(p) =e= Q_sink_d_c(p)/COP_d_c(p)*(COP_d_c(p)+1);
d_mass_flow_rate_source_c(p) .. m_d_source_c(p) =e= Q_source_d_c(p)/cp_d_c(p)/(T_h_o_d_c(p)-T_h_i_d_c(p))*conv;
d_volume_flow_rate_source_c(p,t) .. V_d_source_c(p) =e= m_d_source_c(p)/rho_d_c(p)*3600;

********************************************************************************
* Heat source and heat sink capacity equations
source_cap(p,t) .. Q_source(p,t) =e= V_d_source(p)*rho(p,t)/3600*cp(p,t)*(T_c_i(p,t)-T_c_o(p,t))/conv;
*source_cap(p,t) .. Q_source(p,t) =e= n_bh(p)*V_bh(p)*rho(p,t)/3600*cp(p,t)*(T_c_i(p,t)-T_c_o(p,t))/conv;
supply_cap(p,t)$(COP(p,t)>1) .. Q_sink(p,t) =e= Q_source(p,t)*COP(p,t)/(COP(p,t)-1);
supply_cap_peak('peak',t)$(COP('peak',t)=1) .. Q_sink('peak',t) =e= Q_source('peak',t);

DC_cond_3(p,t) .. Q_free('DC',t) =e= 0 ;
DC_cond_4(p,t) .. Q_hs('DC',t) =e= Q_DC('DC',t)*COP('DC',t)/(COP('DC',t)-1) ;

* names of source and sink were kept from DH, even though meaning is the opposite
* temperatures were switched twice inlet outlet and source sink
source_cap_c(p,t) .. Q_source_c(p,t) =e= V_d_source_c(p)*rho_c(p,t)/3600*cp_c(p,t)*(T_h_o_c(p,t)-T_h_i_c(p,t))/conv;
*source_cap_c(x,t) .. Q_source_c(x,t) =e= n_bh_c(x)*V_bh_c(x)*rho_c(x,t)/3600*cp_c(x,t)*(T_h_o_c(x,t)-T_h_i_c(x,t))/conv;
* changed to +1 for cooling
supply_cap_c(p,t) .. Q_sink_c(p,t) =e= Q_source_c(p,t)*COP_c(p,t)/(COP_c(p,t)+1);

********************************************************************************
*Storage equations
St_char.up(t)= st_char_rate;
St_dis.up(t)= st_char_rate;
Max_storage(t) .. St_lev(t) =l= st_c_max*u_st;
storage_cap(t) .. st_c =g= St_lev(t);
storage_init(t)$(ord(t)=1) .. St_lev(t) =e= st_c_init+St_char(t)-St_dis(t)-St_loss*St_lev(t);
storage(t)$(ord(t)>1) .. St_lev(t) =e= St_lev(t-1)+St_char(t)-St_dis(t)-St_loss*St_lev(t);
storage_end(t)$(ord(t)=%te%) .. St_lev(t) =e= st_c_init;
storage_vol .. st_V =e= st_c/(rho_DH*cp_DH*(T_s_d-T_r_d)/3600/conv);

*Cold Storage equations
St_char_c.up(t)= st_char_rate_c;
St_dis_c.up(t)= st_char_rate_c;
Max_storage_c(t) .. St_lev_c(t) =l= st_c_max_c*u_st_c;
storage_cap_c(t) .. St_c_c =g= St_lev_c(t);
storage_init_c(t)$(ord(t)=1) .. St_lev_c(t) =e= st_c_init_c+St_char_c(t)-St_dis_c(t)-St_loss_c*St_lev_c(t);
storage_c(t)$(ord(t)>1) .. St_lev_c(t) =e= St_lev_c(t-1)+St_char_c(t)-St_dis_c(t)-St_loss_c*St_lev_c(t);
storage_end_c(t)$(ord(t)=%te%) .. St_lev_c(t) =e= st_c_init_c;
storage_vol_c .. st_V_c =e= st_c_c/(rho_DC*cp_DC*(T_r_d_c-T_s_d_c)/3600/conv);

********************************************************************************

check_Q_sink(p,t) .. Q_sink_d_check(p)=g=Q_sink(p,t);
check_Q_hs(p,t) .. Q_hs_check(p)=g=Q_hs(p,t);
check_Q_sink_c(p,t) .. Q_sink_d_check_c(p)=g=Q_sink_c(p,t);
check_Q_hs_c(p,t) .. Q_hs_check_c(p)=g=Q_DC(p,t);

supply(p,t) .. Q_hs(p,t) =l= Q_sink(p,t) ;
supply_max(p,t) .. Q_hs(p,t) =l= Q_sink_d(p) ;
supply_max1(p,t) .. Q_hs(p,t) =l= Q_sink_d_max(p)*u_min(p,t) ;
supply_min(p,t) ..  Q_hs(p,t) =g= 1*u_min(p,t) ;

supply_c(p,t) .. Q_DC(p,t) =l= Q_sink_c(p,t) ;
supply_max_c(p,t) .. Q_DC(p,t) =l= Q_sink_d_c(p) ;
supply_max1_c(p,t) .. Q_DC(p,t) =l= Q_sink_d_max_c(p)*u_min_c(p,t) ;

* free cooling design conditions
Free_power_d(p) .. P_free_d(p) =e= f_aux_el*Q_free_d(p);
Q_free_d.up(p) = 1000;
Q_free_d.up('gw') = 5;
Q_free_d.up('peak') = 0;
Q_sink_d_c.up('peak') = 0;

*free cooling
Scalar T_free_check temnperature at which free cooling will be activated / 3 /;
Scalar T_free_check2 temnperature at which free cooling will be activated / 13 /;
Parameter T_free_o_d(p) design outlet temperature after free cooling;
T_free_o_d(p)=T_free_check2;
Parameter T_free_i_d(p) design inlet source temperature before free cooling;
T_free_i_d(p)=T_free_check;
Parameter T_free_i(p,t) outlet temperature after free cooling;
T_free_i(p,t)=T_h_i_c(p,t);
Parameter T_free_uni(p,t) outlet temperature after free cooling;
T_free_uni(p,t)=T_free_i(p,t);
Parameter T_free_o(p,t) outlet temperature after free cooling;
T_free_o(p,t)$(T_free_i(p,t)<=T_free_check+T_K)=T_h_i_c(p,t)+T_r_c(t)-T_s_c(t);
T_free_o(p,t)$(T_free_i(p,t)>=T_free_check2+T_K)=T_h_i_c(p,t);
T_free_o(p,t)$(T_free_i(p,t)>T_free_check+T_K$(T_free_i(p,t)<T_free_check2+T_K))=T_free_check2+T_K;
Parameter c_free(p) specific investment costs for dry cooler varying with capacity;
c_free('air') = 500*1000/ex_EUR_DKK;
c_free('sea') = 500*1000/ex_EUR_DKK;
c_free('sew') = 500*1000/ex_EUR_DKK;
c_free('gw') = (T_free_check2-T_free_check)/(T_free_check2-T_gw)*500*1000/ex_EUR_DKK;
c_free('DC') = 5000000*1000/ex_EUR_DKK;
c_free('peak') = 5000000*1000/ex_EUR_DKK;

Parameter T_free_o_CK(p,t) outlet temperature after free cooling;
T_free_o_CK(p,t)=T_free_o(p,t)-T_K;

Parameter Q_free_d_max(p) dfe;
Q_free_d_max(p)= 10000;
Free_capacity_d(p,t) .. Q_free(p,t) =l= Q_free_d(p);
Free_capacity_d2(p,t) .. Q_free(p,t) =l= Q_free_d_max(p)*u_free(p) ;
Free_mass_Flow_Rate_max(p,t) .. m_free_d(p) =g= m_free(p,t) ;

Parameter T_test_free(p,t) sd;
T_test_free(p,t)$(T_free_uni(p,t)<T_free_check2+T_K)=T_free_i(p,t);
T_test_free(p,t)$(T_free_uni(p,t)>=T_free_check2+T_K)=0;
Parameter T_s_c_high(p,t) for allowable heat exchange;
T_s_c_high(p,t) = max(T_free_check+T_K+3,min(T_free_i(p,t)+3,T_free_check2+T_K)) ;
Parameter T_s_c_low(p,t) for allowable heat exchange;
T_s_c_low(p,t) = max(T_free_check+T_K,min(T_free_i(p,t),T_free_check2+T_K)) ;

Free_mass_flow_rate_d(p) .. m_free_d(p) =l= Q_free_d(p)/cp_d_c(p)/(T_free_o_d(p)-T_free_i_d(p))*conv;
Free_volume_flow_rate_d(p) .. V_free_d(p) =e=  m_free_d(p)/rho_d_c(p)*3600;

Free_mass_flow_rate(p,t)$(T_free_i(p,t)<=T_free_check+T_K) .. Q_free(p,t) =e= m_free(p,t)*cp_c(p,t)*(T_free_o(p,t)-T_free_i(p,t))/conv;
Free_cooling_x2(p,t)$(T_free_i(p,t)>T_free_check+T_K$(T_free_i(p,t)<T_free_check2+T_K)) .. Q_free(p,t) =l= V_sink_c(p,t)*rho_DC/3600*cp_DC*(T_r_c(t)-T_s_c_high(p,t))/conv;
Free_cooling_x(p,t)$(T_free_i(p,t)>T_free_check+T_K$(T_free_i(p,t)<T_free_check2+T_K)) .. Q_free(p,t) =l= m_free(p,t)*cp_c(p,t)*(T_free_check2+T_K-T_s_c_low(p,t))/conv;

Chiller_cooling(p,t)$(T_free_i(p,t)>=T_free_check2+T_K)  .. m_free(p,t) =e= 0;
Chiller_cooling1(p,t)$(T_free_i(p,t)>=T_free_check2+T_K)  .. Q_free(p,t) =e= 0;

Free_volume_flow_rate(p,t) .. V_free(p,t) =e=  m_free(p,t)/rho_c(p,t)*3600;
Free_power(p,t) .. P_free(p,t) =e= f_aux_el*Q_free(p,t);
Free_power_tot(p) .. P_free_tot(p) =e= sum(t, P_free(p,t)/conv);
Free_cooling_tot(p) .. Q_free_tot(p) =e= sum(t,(Q_free(p,t))/conv);

demand_c(t) .. sum(p, Q_DC(p,t)+Q_free(p,t)) =e= d_DC(t)+St_char_c(t)-St_dis_c(t);
demand(t) .. sum(p, Q_hs(p,t)) =e= d(t)+St_char(t)-St_dis(t)+Q_hs('sew',t)*l_DH;

power_eq(p_nDC,t) .. P_hs(p_nDC,t) =e= Q_hs(p_nDC,t)/COP(p_nDC,t);
power_eq1(p,t) .. P_hs('DC',t) =e= Q_hs('DC',t)/COP('DC',t)/2;
supply_tot(p) .. Q_tot(p) =e= sum(t, Q_hs(p,t))/conv;
power_tot(p) .. P_tot(p) =e= sum(t, P_hs(p,t)/conv);
power_max(p,t) .. P_max(p) =g= P_hs(p,t);

power_eq_c(p_nDC,t) .. P_hs_c(p_nDC,t) =e= Q_DC(p_nDC,t)/COP_c(p_nDC,t);
power_eq_c1(p,t) .. P_hs_c('DC',t) =e= Q_DC('DC',t)/COP_c('DC',t)/2;
supply_tot_c(p) .. Q_tot_c(p) =e= sum(t, Q_DC(p,t))/conv;
power_tot_c(p) .. P_tot_c(p) =e= sum(t, P_hs_c(p,t)/conv);
power_max_c(p,t) .. P_max_c(p) =g= P_hs_c(p,t);

$ontext
*power_aux(p,t) .. P_aux_el_t(p,t) =e= P_aux_el_d(p)*(V_source(p,t)/V_d_source(p))**f_pl;
load_source0(p,t)$(V_d_source(p)>0) .. load_hs(p,t) =e= V_source(p,t)/n_bh(p)/V_bh(p);
load_source(p,t)$(V_d_source(p)>0) .. load_hs(p,t) =e= load_min(p)+sum(k,Pk(p,t,k)) ;
load_source1(p,t)$(V_d_source(p)=0) .. load_hs(p,t) =e= 0;
load_source_max(p,t) .. load_hs_max(p) =g= load_hs(p,t);

Parameter data(k,p,*) p 116;
data(k,p,'DP')$(V_d_source(p)>0) = (load_max(p)-load_min(p))/card(k);
data(k,p,'DP')$(V_d_source(p)=0) =0;
data(k,p,'Pini') = (ord(k)-1)*data(k,p,'DP') + load_min(p);
data(k,p,'Pfin') = data(k,p,'Pini') + data(k,p,'DP');
data(k,p,'Cini') = P_aux_el_d(p)*power(data(k,p,'Pini'),f_pl);
data(k,p,'Cfin') = P_aux_el_d(p)*power(data(k,p,'Pfin'),f_pl);
data(k,p,'s')$(V_d_source(p)>0)    = (data(k,p,'Cfin') - data(k,p,'Cini'))/data(k,p,'DP');
data(k,p,'s')$(V_d_source(p)=0)    = 0;

P_aux_el_t.up(p,t)=P_aux_el_d(p);
P_aux_el_t.lo(p,t)=0;
Pk.up(p,t,k)=data(k,p,'DP');
Pk.lo(p,t,k)=0;

*power_aux(p,t)$(Q_sink_d(p)>0) .. P_aux_el_t(p,t) =e= load_min(p)+sum(k,Pk(p,t,k));
power_aux(p,t)$(Q_sink_d(p)>0) .. P_aux_el_t(p,t) =e= P_aux_el_d(p)*power(load_min(p),f_pl)+sum(k,data(k,p,'s')*Pk(p,t,k));
*load_curve(p,t)$(V_d_source(p)>0) .. load_hs(p,t)=e=
*load_curve1(p,t)$(V_d_source(p)=0) .. load_hs(p,t)=e= 0;
*power_aux(p,t)$(Q_sink_d(p)>0) .. P_aux_el_t(p,t) =e= P_aux_el_d(p)*(load_hs(p,t))**f_pl;


power_aux1(p,t)$(Q_sink_d(p)=0) .. P_aux_el_t(p,t) =e= 0;
power_aux_annual(p) .. power_aux_tot_p(p) =e= sum(t,P_aux_el_t(p,t));

*power_aux(p,t) .. P_aux_el_t(p,t) =e= P_aux_el_d(p)*(V_source(p,t)/V_d_source(p))**f_pl;
load_source0_c(x,t)$(V_d_source_c(x)>0) .. load_hs_c(x,t) =e= V_source_c(x,t)/n_bh_c(x)/V_bh_c(x);
load_source_c(x,t)$(V_d_source_c(x)>0) .. load_hs_c(x,t) =e= load_min_c(x)+sum(k,Pk_c(x,t,k)) ;
load_source1_c(x,t)$(V_d_source_c(x)=0) .. load_hs_c(x,t) =e= 0;
load_source_max_c(x,t) .. load_hs_max_c(x) =g= load_hs_c(x,t);

Parameter data_c(k,x,*) p 116;
data_c(k,x,'DP_c')$(V_d_source_c(x)>0) = (load_max_c(x)-load_min_c(x))/card(k);
data_c(k,x,'DP_c')$(V_d_source_c(x)=0) =0;
data_c(k,x,'Pini_c') = (ord(k)-1)*data_c(k,x,'DP_c') + load_min_c(x);
data_c(k,x,'Pfin_c') = data_c(k,x,'Pini_c') + data_c(k,x,'DP_c');
data_c(k,x,'Cini_c') = P_aux_el_d_c(x)*power(data_c(k,x,'Pini_c'),f_pl);
data_c(k,x,'Cfin_c') = P_aux_el_d_c(x)*power(data_c(k,x,'Pfin_c'),f_pl);
data_c(k,x,'s_c')$(V_d_source_c(x)>0)    = (data_c(k,x,'Cfin_c') - data_c(k,x,'Cini_c'))/data_c(k,x,'DP_c');
data_c(k,x,'s_c')$(V_d_source_c(x)=0)    = 0;

P_aux_el_t_c.up(x,t)=P_aux_el_d_c(x);
P_aux_el_t_c.lo(x,t)=0;
Pk_c.up(x,t,k)=data_c(k,x,'DP_c');
Pk_c.lo(x,t,k)=0;

*power_aux(p,t)$(Q_sink_d(p)>0) .. P_aux_el_t(p,t) =e= load_min(p)+sum(k,Pk(p,t,k));
power_aux_c(x,t)$(Q_sink_d_c(x)>0) .. P_aux_el_t_c(x,t) =e= P_aux_el_d_c(x)*power(load_min_c(x),f_pl)+sum(k,data_c(k,x,'s_c')*Pk_c(x,t,k));
*load_curve(p,t)$(V_d_source(p)>0) .. load_hs(p,t)=e=
*load_curve1(p,t)$(V_d_source(p)=0) .. load_hs(p,t)=e= 0;
*power_aux(p,t)$(Q_sink_d(p)>0) .. P_aux_el_t(p,t) =e= P_aux_el_d(p)*(load_hs(p,t))**f_pl;

power_aux1_c(x,t)$(Q_sink_d_c(x)=0) .. P_aux_el_t_c(x,t) =e= 0;
power_aux_annual_c(x) .. power_aux_tot_p_c(x) =e= sum(t,P_aux_el_t_c(x,t));
$offtext

********************************************************************************
*volume flow rates
real_source_cap(p,t) .. Q_source_real(p,t) =e= Q_hs(p,t)/COP(p,t)*(COP(p,t)-1);
volume_flow_rate_source(p,t) .. V_source(p,t) =e= Q_source_real(p,t)/(rho(p,t)/3600*cp(p,t)*(T_c_i(p,t)-T_c_o(p,t))/conv);
V_source.up(p,t) = V_d_source_max(p,t);
volume_flow_rate_source_max(p,t) .. V_source(p,t) =l= V_d_source(p);

volume_flow_rate_sink(p,t) .. V_sink(p,t) =e= Q_hs(p,t)/(rho_DH/3600*cp_DH*(T_s(t)-T_r(t))/conv);
volume_flow_rate_sink_max(p,t) .. V_sink(p,t) =l= V_DH_max;
V_sink.up(p,t) = V_DH_max;
volume_flow_rate_sink_max_real(p,t) .. V_sink_max(p) =g= V_sink(p,t);
*Diameter_source(p).. d_DH(p) =e= sqrt(V_sink_max(p)/(u_DH_max*pi/4*3600));

*volume flow rates
real_source_cap_c(p,t) .. Q_source_real_c(p,t) =e= Q_DC(p,t)/COP_c(p,t)*(COP_c(p,t)+1)+Q_free(p,t);
volume_flow_rate_source_c(p,t) .. V_source_c(p,t) =e= (Q_source_real_c(p,t))/(rho_c(p,t)/3600*cp_c(p,t)*(T_h_o_c(p,t)-T_h_i_c(p,t))/conv);
V_source_c.up(p,t) = V_d_source_max_c(p,t);
volume_flow_rate_source_max_c(p,t) .. V_source_c(p,t) =l= V_d_source_max_c(p,t);

volume_flow_rate_sink_c(p,t) .. V_sink_c(p,t) =e= (Q_DC(p,t)+Q_free(p,t))/(rho_DC/3600*cp_DC*(T_r_c(t)-T_s_c(t))/conv);
volume_flow_rate_sink_max_c(p,t) .. V_sink_c(p,t) =l= V_DH_max_c;
V_sink_c.up(p,t) = V_DH_max_c;
volume_flow_rate_sink_max_real_c(p,t) .. V_sink_max_c(p) =g= V_sink_c(p,t);

********************************************************************************
* investment costs heat sources
*inv_tot_base(p) .. C_inv_tot('base') =e= c_inv_min('base')+ c_inv_f('base')*Q_sink_d('base');
*inv_tot_peak(p) .. C_inv_tot('peak') =e= c_inv_min('peak')+ c_inv_f('peak')*Q_sink_d('peak');
inv_tot_hs(p) .. C_inv_tot(p) =e= c_inv_min0(p)*u(p)+ c_inv_f0(p)*Q_sink_d(p);
inv_tot_hs_c(p) .. C_inv_tot_c(p) =e= c_inv_min_c0(p)*u_c(p) + c_inv_f_c(p)*Q_sink_d_c(p);

* investment costs storage
inv_stor .. C_inv_stor =e= (c_inv_stor_f*st_V+u_st*c_inv_stor_min)*a_c_inv_st/a_c_inv;
inv_stor_c .. C_inv_stor_c =e= (c_inv_stor_f_c*st_V_c)*a_c_inv_st/a_c_inv;

* for 2.0 m/s
inv_DH(pns) .. C_DH(pns) =e= ((3.1*V_sink_max(pns)+505*u(pns))*pipe_l(pns))*a_c_inv_DH/a_c_inv;

* double check the numbers and check here depending on the purpose whether investment counts once or twice
Parameter u_pipe(p) sd;
u_pipe(p)=1;
u_pipe('air')=0;
u_pipe('gw')=0;
u_pipe('DC')=0;

* check whether investment is accounted for twice
inv_DH_c(pns) .. C_DH_c(pns) =e= (500/ex_EUR_DKK*conv*Q_sink_d_c(pns)*u_pipe(pns)+3000/ex_EUR_DKK*u_c(pns)*pipe_l_c(pns))*a_c_inv_DH/a_c_inv;
inv_DH_gw_c(p) .. C_DH_c('gw') =e= (500/ex_EUR_DKK*conv*Q_sink_d_c('gw')*u_pipe('gw')+3000/ex_EUR_DKK*u_c('gw')*pipe_l_c('gw'))*a_c_inv_DH/a_c_inv;
inv_DH_sea_c(p) .. C_DH_c('sea') =e= (500/ex_EUR_DKK*conv*Q_sink_d_c('sea')*u_pipe('sea')+3000/ex_EUR_DKK*u_c('sea')*pipe_l_c('sea'))*a_c_inv_DH/a_c_inv;

inv_DH_sea(p) .. C_DH('sea') =e= (161100*Q_sink_d('sea')+u('sea')*3758400)*a_c_inv_DH/a_c_inv;
*inv_DH_gw(p) .. C_DH('gw') =e= ((3.1*V_d_source('gw')+505*u('gw'))*pipe_l('gw'))*a_c_inv_DH/a_c_inv;
inv_DH_gw(p) .. C_DH('gw') =e= (500/ex_EUR_DKK*conv*Q_sink_d('gw')*u_pipe('gw')+3000/ex_EUR_DKK*u('gw')*pipe_l('gw'))*a_c_inv_DH/a_c_inv;

* check whether investment is accounted for twice
*free cooling
u_z_free.lo(p)=0;
u_z_free.up(p)=1;
u_z_free_Eq1(p) .. u_z_free(p) =l= u_free(p);
u_z_free_Eq2(p) .. u_z_free(p) =l= 1-u_c(p);
u_z_free_Eq3(p) .. u_z_free(p) =g= u_free(p)-u_c(p);
*u_z_free(p)
inv_free2(p) .. C_inv_free(p) =e= c_free(p)*Q_free_d(p)+c_inv_min_free(p)*u_free(p)+u_z_free(p)*3000/ex_EUR_DKK*pipe_l_c(p)*a_c_inv_DH/a_c_inv;
inv_free_tot .. C_inv_free_tot =e=  sum(p,C_inv_free(p));


********************************************************************************
* Maintenance costs
*maintenance_base(p) .. C_m('base') =e= m_f('base')*Q_tot('base')+m_min('base')*Q_sink_d('base');
*maintenance_peak(p) .. C_m('peak') =e= m_f('peak')*Q_tot('peak')+m_min('peak')*Q_sink_d('peak');
maintenance_hs(p_nDC) .. C_m(p_nDC) =e= m_f(p_nDC)*Q_tot(p_nDC)+m_min(p_nDC)*Q_sink_d(p_nDC);
maintenance_hs1(p) .. C_m('DC') =e= (m_f('DC')*Q_tot('DC')+m_min('DC')*Q_sink_d('DC'))/2;

maintenance_hs_c(p_nDC) .. C_m_c(p_nDC) =e=  m_min_c(p_nDC)*(C_inv_tot_c(p_nDC)+C_inv_free(p_nDC));
maintenance_hs_c1(p) .. C_m_c('DC') =e= (m_f('DC')*Q_tot('DC')+m_min('DC')*Q_sink_d('DC'))/2;

********************************************************************************

*cost function for each production unit
annual_costs_el_prod(p) .. C_el_p(p) =e= sum((t),P_hs(p,t)*(p_el(t)+tax_c_all(t)));
Total_annual_costs_el_prod .. C_el_tot =e= sum(p,C_el_p(p));
Total_costs_investments .. C_inv_tot_all =e= sum(p,C_inv_tot(p));
Total_annual_costs_maintenance .. C_m_tot =e= sum(p,C_m(p));
Total_investment_costs_DH_pipes .. C_DH_pipe =e= sum(p,C_DH(p));

annual_costs_el_prod_c(p) .. C_el_p_c(p) =e= sum((t),(P_hs_c(p,t)+P_free(p,t))*(p_el(t)+tax_c_all(t)));
Total_annual_costs_el_prod_c .. C_el_tot_c =e= sum(p,C_el_p_c(p));
Total_costs_investments_c .. C_inv_tot_all_c =e= sum(p,C_inv_tot_c(p));
Total_annual_costs_maintenance_c .. C_m_tot_c =e= sum(p,C_m_c(p));
Total_investment_costs_DH_pipes_c .. C_DH_pipe_c =e= sum(p,C_DH_c(p));

*aux el costs
**Aux_el_costs(p) .. C_el_p_aux(p) =e= sum(t,(P_aux_el_t(p,t)*(p_el(t)+(fee_trans+fee_distr(t))/ex_EUR_DKK)));
**Aux_el_costs1 .. C_el_aux =e= sum(p,C_el_p_aux(p));
**Aux_el_costs_c(x) .. C_el_p_aux_c(x) =e= sum(t,(P_aux_el_t_c(x,t)*(p_el(t)+(fee_trans+fee_distr(t))/ex_EUR_DKK)));
**Aux_el_costs1_c .. C_el_aux_c =e= sum(x,C_el_p_aux_c(x));

*cost function
********************************************************************************
cost_p(p) .. Z_p(p) =e= C_el_p(p) + (C_inv_tot(p) + C_DH(p))*a_c_inv + C_m(p);
cost_p_c(p) .. Z_p_c(p) =e= C_el_p_c(p) + (C_inv_tot_c(p) + C_DH_c(p))*a_c_inv + C_m_c(p)+C_inv_free(p)*a_c_inv ;
**cost_p(p) .. Z_p(p) =e= C_el_p(p) + (C_inv_tot(p) + C_DH(p))*a_c_inv + C_m(p) + C_el_p_aux(p);
**cost_p_c(x) .. Z_p_c(x) =e= C_el_p_c(x) + (C_inv_tot_c(x) + C_DH_c(x))*a_c_inv + C_m_c(x) + C_el_p_aux_c(x);

*cost .. Z =e= (sum(p,(C_el_p(p) + (C_inv_tot(p) + C_DH(p))*a_c_inv + C_m(p)))+C_inv_stor*a_c_inv+C_inv_DH_net*a_c_inv);
cost .. Z =e= (sum(p,(C_el_p(p) + (C_inv_tot(p) + C_DH(p))*a_c_inv + C_m(p)))+C_inv_stor*a_c_inv+C_inv_DH_net*a_c_inv)+
(sum(p,(C_el_p_c(p) + (C_inv_tot_c(p) + C_DH_c(p))*a_c_inv + C_m_c(p)))+C_inv_stor_c*a_c_inv+C_inv_free_tot*a_c_inv);


* minimize CO2 emissions
* CO2_HP_c = sum((p,t),(P_hs_c.l(p,t)+P_free.l(p,t))*CO2(t))/sum(t,d_DC(t));
* CO2_HP = sum((p,t),(P_hs.l(p,t))*CO2(t))/sum(t,d(t));
Parameter d_tot total heat demand;
d_tot=sum(t,d(t));
Parameter d_tot_C total heat demand;
d_tot_c=sum(t,d_DC(t));

Variables
Z_CO2
Z_CO2_t(p,t)
Z_CO2_p(p)
;

*Positive Variable Z_CO2, Z_CO2_p,Z_CO2_t;

Equations
emissions
emissions_t(p,t)
emissions_p(p)
;
emissions_t(p,t) .. Z_CO2_t(p,t) =e= CO2(t)*P_hs_c(p,t)+CO2(t)*P_free(p,t)+CO2(t)*P_hs(p,t);
emissions_p(p) .. Z_CO2_p(p) =e= sum((t),Z_CO2_t(p,t))/conv;
emissions .. Z_CO2 =e= sum((p),Z_CO2_p(p));
*Z_CO2.up=4840;

* maximize performance
Variables
Z_COP
Z_COP_p(p)
*P_tot_all
;

*Positive Variable Z_COP, Z_COP_p;

Equations
performance
performance_t(p)
;
performance_t(p) .. Z_COP_p(p) =e= sum(t,P_hs_c(p,t)+P_free(p,t)+P_hs(p,t))/conv;
*performance_t(p) .. Z_COP_p(p) =e= 0;
performance .. Z_COP =e= sum((p),Z_COP_p(p));
*Z_COP.up=25.870;

**cost .. Z =e= (sum(p,(C_el_p(p) + (C_inv_tot(p) + C_DH(p))*a_c_inv + C_m(p)) + C_el_p_aux(p))+C_inv_stor*a_c_inv+C_inv_DH_net*a_c_inv)+
**(sum(x,(C_el_p_c(x) + (C_inv_tot_c(x) + C_DH_c(x))*a_c_inv + C_m_c(x)) + C_el_p_aux_c(x))+C_inv_stor_c*a_c_inv);

********************************************************************************

Parameter Q_sink_d_max(p) ;
Q_sink_d_max('DC')=0;

*relative gap
option optcr = 0.001
*option optcr = 0.00001
option resLim = 12000


Model DHC_separate / all /;

Solve DHC_separate using mip minimizing Z ;
*Solve DHC_separate using mip minimizing Z_CO2 ;

*Solve DHC_separate using mip minimizing Z_CO2 ;
*Solve DHC_separate using mip maximizing Z_COP ;
*Solve heatsources using miqcp minimizing Z ;

********************************************************************************
********************************************************************************

Parameter Q_sink_d_max(p) ;
Q_sink_d_max('DC')=1000;


***************************************************************
Variables

u_z(p) multiply binary variables u(p) and u_c(p)
u_z_c(p) multiply binary variables u(p) and u_c(p)
;

Positive Variable u_z, u_z_c;

Equations

u_z_Eq1(p)
u_z_Eq2(p)
u_z_Eq3(p)
u_z_c_Eq1(p)
u_z_c_Eq2(p)
u_z_c_Eq3(p)
inv_tot_hs0(p)
inv_tot_hs_c0(p_nDC)
inv_tot_hs_c01(p)
inv_DH_gw_c0(p)
inv_DH_sea_c0(p)
;

u_z.lo(p)=0;
u_z.up(p)=1;
u_z_Eq1(p) .. u_z(p) =l= u_c(p);
u_z_Eq2(p) .. u_z(p) =l= 1-u(p);
u_z_Eq3(p) .. u_z(p) =g= u_c(p)-u(p);

u_z_c.lo(p)=0;
u_z_c.up(p)=1;
u_z_c_Eq1(p) .. u_z_c(p) =l= u(p);
u_z_c_Eq2(p) .. u_z_c(p) =l= 1-u_c(p);
u_z_c_Eq3(p) .. u_z_c(p) =g= u(p)-u_c(p);

*inv_tot_hs0(p) .. C_inv_tot(p) =e= c_inv_min(p)*u(p)+c_inv_min2(p)*u_z_c(p)+ c_inv_f(p)*Q_sink_d(p)+c_inv_f2(p)*Q_DH_req(p);
inv_tot_hs0(p) .. C_inv_tot(p) =e= c_inv_min(p)*u(p)+c_inv_min2(p)*u_z_c(p)+ c_inv_f0(p)*Q_sink_d(p);

inv_tot_hs_c0(p_nDC) .. C_inv_tot_c(p_nDC) =e= u_z(p_nDC)*c_inv_min_c(p_nDC)+u_c(p_nDC)*c_inv_min_c2(p_nDC) + c_inv_f_c(p_nDC)*Q_sink_d_c(p_nDC);
inv_tot_hs_c01(p) .. C_inv_tot_c('DC') =e= u_c('DC')*c_inv_min_c('DC')+u_c('DC')*c_inv_min_c2('DC') + c_inv_f_c('DC')*Q_sink_d_c('DC');

inv_DH_gw_c0(p) .. C_DH_c('gw') =e= (500/ex_EUR_DKK*conv*Q_sink_d_c('gw')*u_pipe('gw')+3000/ex_EUR_DKK*u_z('gw')*pipe_l_c('gw'))*a_c_inv_DH/a_c_inv;
inv_DH_sea_c0(p) .. C_DH_c('sea') =e= (500/ex_EUR_DKK*conv*Q_sink_d_c('sea')*u_pipe('sea')+3000/ex_EUR_DKK*u_z('sea')*pipe_l_c('sea'))*a_c_inv_DH/a_c_inv;

*cost .. Z =e= (sum(p,(C_el_p(p) + (C_inv_tot(p) + C_DH(p))*a_c_inv + C_m(p)))+C_inv_stor*a_c_inv+C_inv_DH_net*a_c_inv)+
*(sum(p,(C_el_p_c(p) + (C_inv_tot_c(p) + C_DH_c(p))*a_c_inv + C_m_c(p)))+C_inv_stor_c*a_c_inv+C_inv_free_tot*a_c_inv);

*relative gap
option optcr = 0.003
*option optcr = 0.00001
option resLim = 18000

set counter / c1*c7 /;
Scalar E;
Parameter report(counter,*);
Parameter report_p(counter,*,p);

Model DHC_with_shared_inv / all - inv_tot_hs-inv_tot_hs_c-inv_tot_hs_c-inv_DH_gw_c-inv_DH_sea_c /;

Parameter ranges(*);
Solve DHC_with_shared_inv using mip minimizing Z ;
*Solve DHC_with_shared_inv using mip minimizing Z_CO2 ;
*Solve DHC_with_shared_inv using mip maximizing Z_COP ;

Parameter f_hs(p) heat sources shares between HPs;
f_hs(p)$(Q_sink_d.l(p)>0) = Q_sink_d.l(p)/(Q_sink_d.l('air')+Q_sink_d.l('gw')+Q_sink_d.l('sea')+Q_sink_d.l('sew')+Q_sink_d.l('peak')+Q_sink_d.l('DC'));
f_hs(p)$(Q_sink_d.l(p)=0)= 0;
Parameter SCOP(p) seasonal COP of each production unit used to cover demand;
SCOP(p_nDC)$(P_tot.l(p_nDC)>0)= Q_tot.l(p_nDC)/(P_tot.l(p_nDC));
SCOP('DC')$(P_tot.l('DC')>0)= Q_tot.l('DC')/(P_tot.l('DC')*2);
*SCOP(p)$(P_tot.l(p)>0)= Q_tot.l(p)/(P_tot.l(p) + Power_aux_tot_p.l(p)/conv);
SCOP(p)$(P_tot.l(p)=0)= 0;
Parameter SCOP_tot seasonal COP of all production units used together to cover demand;
*SCOP_tot= sum(p,Q_tot.l(p))/sum(p,P_tot.l(p));
SCOP_tot = sum(t,d(t))/conv/(sum(p_nDC,P_tot.l(p_nDC))+P_tot.l('DC')*2);

Parameter C_prod_HP_inv Costs of HP production;
*C_prod_HP_inv = ((C_inv_tot_all.l)*a_c_inv )/(sum(p,Q_tot.l(p))*conv);
C_prod_HP_inv = ((C_inv_tot_all.l)*a_c_inv )/sum(t,d(t));
* storage costs
Parameter C_inv_stor_p(p) total investment costs for storage for each heat source;
*C_inv_stor_p(p) = C_inv_stor.l*Q_sink_d(p)/(Q_sink_d('air')+Q_sink_d('gw')+Q_sink_d('sea')+Q_sink_d('sew')+Q_sink_d('peak'));
C_inv_stor_p(p)$(Q_tot.l(p)>0) = (C_inv_stor.l*a_c_inv)*Q_sink_d.l(p)/(Q_sink_d.l('air')+Q_sink_d.l('gw')+Q_sink_d.l('sea')+Q_sink_d.l('sew')+Q_sink_d.l('peak')+Q_sink_d.l('DC'))/(Q_tot.l(p)*conv);
C_inv_stor_p(p)$(Q_tot.l(p)=0) = 0;
*C_inv_stor_p(p) = C_inv_stor.l*Q_sink_d(p)/sum(p,Q_sink_d(p));
Parameter C_prod_HP_st Costs of HP production;
*C_prod_HP_st = ((C_inv_stor.l)*a_c_inv)/(sum(p,Q_tot.l(p))*conv);
C_prod_HP_st = ((C_inv_stor.l)*a_c_inv)/sum(t,d(t));
* DH netowork costs
Parameter C_prod_DH_net_p(p) total investment costs for storage for each heat source;
C_prod_DH_net_p(p)$(Q_tot.l(p)>0) = (C_inv_DH_net*a_c_inv)*Q_sink_d.l(p)/(Q_sink_d.l('air')+Q_sink_d.l('gw')+Q_sink_d.l('sea')+Q_sink_d.l('sew')+Q_sink_d.l('peak')+Q_sink_d.l('DC'))/(Q_tot.l(p)*conv);
C_prod_DH_net_p(p)$(Q_tot.l(p)=0) = 0;
Parameter C_prod_DH_net Costs of DH network;
C_prod_DH_net = ((C_inv_DH_net)*a_c_inv)/sum(t,d(t));

* Heat source investment costs
Parameter C_prod_HP_DH_p(p) Costs of HP production;
C_prod_HP_DH_p(p)$(Q_tot.l(p)>0) = (C_DH.l(p)*a_c_inv)/(Q_tot.l(p)*conv);
C_prod_HP_DH_p(p)$(Q_tot.l(p)=0) = 0;
Parameter C_prod_HP_DH Costs of HP production;
*C_prod_HP_DH = ((C_DH_pipe.l)*a_c_inv)/(sum(p,Q_tot.l(p))*conv);
C_prod_HP_DH = ((C_DH_pipe.l)*a_c_inv)/sum(t,d(t));
* Maintenance costs
Parameter C_prod_HP_m_p(p) Costs of HP production;
C_prod_HP_m_p(p)$(Q_tot.l(p)>0) = (C_m.l(p))/(Q_tot.l(p)*conv);
C_prod_HP_m_p(p)$(Q_tot.l(p)=0) = 0;
Parameter C_prod_HP_m Costs of HP production;
*C_prod_HP_m = (C_m_tot.l)/(sum(p,Q_tot.l(p))*conv);
C_prod_HP_m = (C_m_tot.l)/sum(t,d(t));
* El costs
Parameter C_prod_HP_el_p(p) Costs of HP production;
C_prod_HP_el_p(p)$(Q_tot.l(p)>0) = (C_el_p.l(p))/(Q_tot.l(p)*conv);
C_prod_HP_el_p(p)$(Q_tot.l(p)=0) = 0;
Parameter C_prod_HP_el Costs of HP production;
C_prod_HP_el = (C_el_tot.l)/(sum(t,d(t)));
* Cost of HP production
Parameter C_prod_HP_p(p) Costs of HP production;
C_prod_HP_p(p)$(Q_tot.l(p)>0) = ((C_inv_tot.l(p)+C_DH.l(p))*a_c_inv  + C_m.l(p) + C_el_p.l(p) )/(Q_tot.l(p)*conv)+C_inv_stor_p(p)+C_prod_DH_net_p(p);
*C_prod_HP_p(p)$(Q_tot.l(p)>0) = ((C_inv_tot.l(p)+C_DH.l(p))*a_c_inv  + C_m.l(p) + C_el_p.l(p) + C_el_p_aux.l(p))/(Q_tot.l(p)*conv)+C_inv_stor_p(p)+C_prod_DH_net_p(p);
C_prod_HP_p(p)$(Q_tot.l(p)=0) = 0;
Parameter C_prod_HP Costs of HP production;
*C_prod_HP = ((C_inv_tot_all.l+C_inv_stor.l+C_DH_pipe.l)*a_c_inv + C_m_tot.l + C_el_tot.l)/(sum(p,Q_tot.l(p))*conv);
C_prod_HP = ((C_inv_tot_all.l+C_inv_stor.l+C_inv_DH_net+C_DH_pipe.l)*a_c_inv + C_m_tot.l + C_el_tot.l)/sum(t,d(t));
*C_prod_HP = ((C_inv_tot_all.l+C_inv_stor.l+C_inv_DH_net+C_DH_pipe.l)*a_c_inv + C_m_tot.l + C_el_tot.l + C_el_aux.l)/sum(t,d(t));

Parameter CO2_HP_p(p) CO2 emissions for using HPs to supply heat;
CO2_HP_p(p)$(Q_tot.l(p)>0) = sum(t,(P_hs.l(p,t))*CO2(t)/conv)/Q_tot.l(p) ;
CO2_HP_p(p)$(Q_tot.l(p)=0)= 0;
*CO2_HP_p(p) = sum(t,(P_hs.l(p,t)+P_aux_el_t.l(p,t))*CO2(t)/conv) ;
Parameter CO2_HP CO2 emissions for using HPs to supply heat;
*CO2_HP = sum((p,t),P_hs.l(p,t)*CO2(t))/sum(p,(Q_tot.l(p)*conv)) ;
CO2_HP = sum((p,t),(P_hs.l(p,t))*CO2(t))/sum(t,d(t));

Parameter f_DC(p) heat sources shares between HPs;
f_DC(p)$(Q_sink_d_c.l(p)>0) = (Q_sink_d_c.l(p)+Q_free_d.l(p))/(Q_sink_d_c.l('air')+Q_free_d.l('air')+Q_sink_d_c.l('sea')+Q_free_d.l('sea')+Q_sink_d_c.l('sew')+Q_free_d.l('sew')+Q_sink_d_c.l('gw')+Q_free_d.l('gw')+Q_free_d.l('DC'));
f_DC(p)$(Q_sink_d_c.l(p)=0)= 0;
Parameter SCOP_c(p) seasonal COP of each production unit used to cover demand;
SCOP_c(p_nDC)$(P_tot_c.l(p_nDC)>0)= (Q_tot_c.l(p_nDC)+Q_free_tot.l(p_nDC))/(P_tot_c.l(p_nDC)+P_free_tot.l(p_nDC) );
SCOP_c('DC')$(P_tot_c.l('DC')>0)= (Q_tot_c.l('DC')+Q_free_tot.l('DC'))/(P_tot_c.l('DC')*2+P_free_tot.l('DC') );
*SCOP_c(x)$(P_tot_c.l(x)>0)= Q_tot_c.l(x)/(P_tot_c.l(x) + Power_aux_tot_p_c.l(x)/conv);
SCOP_c(p)$(P_tot_c.l(p)=0)= 0;
Parameter SCOP_tot_c seasonal COP of all production units used together to cover demand;
*SCOP_tot= sum(p,Q_tot.l(p))/sum(p,P_tot.l(p));
SCOP_tot_c = sum(t,d_DC(t))/conv/(sum(p_nDC,P_tot_c.l(p_nDC)+P_free_tot.l(p_nDC))+P_tot_c.l('DC')*2+P_free_tot.l('DC'));

Parameter SCOP_sys seasonal COP considering heating and cooling;
SCOP_sys=  sum(t,d_DC(t)+d(t))/conv/sum(p,P_tot_c.l(p)+P_free_tot.l(p)+P_tot.l(p));

* HP investment costs

Parameter C_prod_HP_inv_p_c(p) Costs of HP production;
C_prod_HP_inv_p_c(p)$(Q_tot_c.l(p)>0) = ((C_inv_tot_c.l(p)+C_inv_free.l(p))*a_c_inv)/((Q_tot_c.l(p)+Q_free_tot.l(p))*conv);
C_prod_HP_inv_p_c(p)$(Q_tot_c.l(p)=0) = 0;
Parameter C_prod_HP_inv_c Costs of HP production;
C_prod_HP_inv_c = ((C_inv_tot_all_c.l+C_inv_free_tot.l)*a_c_inv )/sum(t,d_DC(t));
* storage costs
Parameter C_inv_stor_p_c(p) total investment costs for storage for each heat source;
* requires update when more heat sources added
C_inv_stor_p_c(p)$(Q_tot_c.l(p)>0) = (C_inv_stor_c.l*a_c_inv)*(Q_sink_d_c.l(p)+Q_free_d.l(p))/(Q_sink_d_c.l('air')+Q_sink_d_c.l('gw')+Q_sink_d_c.l('sea')+Q_sink_d_c.l('sew')+Q_sink_d_c.l('DC')+Q_free_d.l('air')+Q_free_d.l('gw')+Q_free_d.l('sea')+Q_free_d.l('sew')+Q_free_d.l('DC'))/((Q_tot_c.l(p)+Q_free_tot.l(p))*conv);
C_inv_stor_p_c(p)$(Q_tot_c.l(p)=0) = 0;
Parameter C_prod_HP_st_c Costs of HP production;
C_prod_HP_st_c = ((C_inv_stor_c.l)*a_c_inv)/sum(t,d_DC(t));
* Heat source investment costs
Parameter C_prod_HP_DH_p_c(p) Costs of HP production;
C_prod_HP_DH_p_c(p)$(Q_tot_c.l(p)>0) = (C_DH_c.l(p)*a_c_inv)/((Q_tot_c.l(p)+Q_free_tot.l(p))*conv);
C_prod_HP_DH_p_c(p)$(Q_tot_c.l(p)=0) = 0;
Parameter C_prod_HP_DH_c Costs of HP production;
C_prod_HP_DH_c = ((C_DH_pipe_c.l)*a_c_inv)/sum(t,d_DC(t));
* Maintenance costs
Parameter C_prod_HP_m_p_c(p) Costs of HP production;
C_prod_HP_m_p_c(p)$(Q_tot_c.l(p)>0) = (C_m_c.l(p))/((Q_tot_c.l(p)+Q_free_tot.l(p))*conv);
C_prod_HP_m_p_c(p)$(Q_tot_c.l(p)=0) = 0;
Parameter C_prod_HP_m_c Costs of HP production;
C_prod_HP_m_c = (C_m_tot_c.l)/sum(t,d_DC(t));
* Aux el consumption costs
*Parameter C_prod_HP_aux_p_c(x) Costs of HP production;
*C_prod_HP_aux_p_c(x)$(Q_tot_c.l(x)>0) = (C_el_p_aux_c.l(x))/(Q_tot_c.l(x)*conv);
*C_prod_HP_aux_p_c(x)$(Q_tot_c.l(x)=0) = 0;
*Parameter C_prod_HP_aux_c Costs of HP production;
*C_prod_HP_aux_c = (C_el_aux_c.l)/(sum(t,d_DC(t)));
* El costs
Parameter C_prod_HP_el_p_c(p) Costs of HP production;
C_prod_HP_el_p_c(p)$(Q_tot_c.l(p)>0) = (C_el_p_c.l(p))/((Q_tot_c.l(p)+Q_free_tot.l(p))*conv);
C_prod_HP_el_p_c(p)$(Q_tot_c.l(p)=0) = 0;
Parameter C_prod_HP_el_c Costs of HP production;
C_prod_HP_el_c = (C_el_tot_c.l)/(sum(t,d_DC(t)));

Parameter C_prod_free_p_c(p) Costs of HP production;
C_prod_free_p_c(p)$(Q_free_tot.l(p)>0) = (C_inv_free.l(p)*a_c_inv)/((Q_tot_c.l(p)+Q_free_tot.l(p))*conv);
C_prod_free_p_c(p)$(Q_free_tot.l(p)=0) = 0;
Parameter C_prod_free_c Costs of HP production;
C_prod_free_c = (C_inv_free_tot.l*a_c_inv)/(sum(t,d_DC(t)));
*inv_free_tot .. C_inv_free_tot =e=  sum(p,C_inv_free(p));

* Cost of HP production
Parameter C_prod_HP_p_c(p) Costs of HP production;
C_prod_HP_p_c(p)$(Q_tot_c.l(p)>0) = ((C_inv_tot_c.l(p)+C_DH_c.l(p)+C_inv_free.l(p))*a_c_inv  + C_m_c.l(p) + C_el_p_c.l(p) )/((Q_tot_c.l(p)+Q_free_tot.l(p))*conv)+C_inv_stor_p_c(p);
*C_prod_HP_p_c(x)$(Q_tot_c.l(x)>0) = ((C_inv_tot_c.l(x)+C_DH_c.l(x))*a_c_inv  + C_m_c.l(x) + C_el_p_c.l(x) + C_el_p_aux_c.l(x))/(Q_tot_c.l(x)*conv)+C_inv_stor_p_c(x);
C_prod_HP_p_c(p)$(Q_tot_c.l(p)=0) = 0;
Parameter C_prod_HP_c Costs of HP production;
C_prod_HP_c = ((C_inv_tot_all_c.l+C_inv_stor_c.l+C_DH_pipe_c.l+C_inv_free_tot.l)*a_c_inv + C_m_tot_c.l + C_el_tot_c.l )/sum(t,d_DC(t));

Parameter CO2_HP_p_c(p) CO2 emissions for using HPs to supply heat;
CO2_HP_p_c(p)$(Q_tot_c.l(p)>0) = sum(t,(P_hs_c.l(p,t)+P_free.l(p,t))*CO2(t)/conv)/(Q_tot_c.l(p)+Q_free_tot.l(p)) ;
CO2_HP_p_c(p)$(Q_tot_c.l(p)=0)= 0;
Parameter CO2_HP_c CO2 emissions for using HPs to supply heat;
*CO2_HP = sum((p,t),P_hs.l(p,t)*CO2(t))/sum(p,(Q_tot.l(p)*conv)) ;
CO2_HP_c = sum((p,t),(P_hs_c.l(p,t)+P_free.l(p,t))*CO2(t))/sum(t,d_DC(t));

*Economic results
Display Z.l, Z_CO2.l,Z_COP.l,C_prod_HP,C_prod_HP_inv,C_prod_HP_st,C_prod_HP_DH,C_prod_HP_m,C_prod_HP_el;
Display Z_p.l,Z_CO2_p.l,Z_COP_p.l, C_prod_HP_p,C_inv_stor_p,C_prod_DH_net_p,C_prod_HP_DH_p,C_prod_HP_m_p;
*Technical results
Display Q_tot.l, P_tot.l, f_hs,d_tot, Q_sink_d.l,SCOP, SCOP_tot, COP_d, eta_L_d;
*Environmental results
Display CO2_HP, CO2_HP_p;
*Remaining results
Display Q_hs.l, P_hs.l, V_d_source.l, V_source.l, V_DH_max, V_sink_max.l, St_V.l, St_c.l, St_lev.l, St_char.l, St_dis.l;
*Economic results
Display C_prod_HP_c,C_prod_HP_inv_c,C_prod_free_c,C_prod_HP_st_c,C_prod_HP_DH_c,C_prod_HP_m_c,C_prod_HP_el_c ;
Display Z_p_c.l,C_prod_HP_p_c,C_prod_HP_inv_p_c,C_prod_free_p_c,C_inv_stor_p_c,C_prod_HP_DH_p_c,C_prod_HP_m_p_c;
*Technical results
Display f_DC,d_tot_c, Q_sink_d_c.l, Q_free_d.l,SCOP_c, SCOP_tot_c, SCOP_sys,COP_d_c, eta_L_d_c;
Display Q_tot_c.l, P_tot_c.l, Q_sink_d_check_c.l, P_free_tot.l, Q_free_tot.l;
*Environmental results
Display CO2_HP_c, CO2_HP_p_c;
*Remaining results
Display Q_DC.l,Q_free.l, P_hs_c.l,V_source_c.l,V_sink_max_c.l,St_V_c.l, St_c_c.l, St_lev_c.l, St_char_c.l, St_dis_c.l;
Display V_d_source_c.l,V_free_d.l, V_DH_max_c, V_sink_c.l, d_DC,C_inv_tot.l,C_inv_tot_c.l, u_z_c.l, u_z.l, u.l, u_c.l;


Parameter SCOP_coef;
SCOP_coef=(sum(p_nDC,P_tot.l(p_nDC))+P_tot.l('DC')*2);
Parameter SCOP_coef1;
SCOP_coef1=(sum(p_nDC,P_tot_c.l(p_nDC)+P_free_tot.l(p_nDC))+P_tot_c.l('DC')*2+P_free_tot.l('DC'));
Parameter SCOP_coef2;
SCOP_coef2=sum(p,P_tot_c.l(p)+P_free_tot.l(p)+P_tot.l(p));

EXECUTE_UNLOAD 'GAMS_output_DHC_with_shared_inv', CO2_HP,f_hs, St_c, St_V, P_hs, Q_hs, P_tot, Q_tot, Q_source, Q_sink, St_lev, St_char, St_dis, COP, f_hs, T_r, T_s, Z,COP, eta_L, C_prod_HP,C_prod_HP_inv, C_prod_HP_st, C_prod_HP_DH, C_prod_HP_m, SCOP_tot,Q_DC,Q_free,St_lev_c, St_char_c, St_dis_c;


*GAMS ELN_opt_DHC_all_free_DCc_CO2 save=DHC_with_shared_inv_results
$ontext

ranges('OF1min')=Z.l;
ranges('OF2max')=Z_CO2.l;
*ranges('OF2max')=Z_COP.l;

Solve DHC_with_shared_inv using mip minimizing Z_CO2 ;
Display Z.l,Z_CO2.l;
*Solve DHC_with_shared_inv using mip minimizing Z_COP ;
*Display Z.l,Z_COP.l;

*$ontext
ranges('OF2min')=Z_CO2.l;
*ranges('OF2min')=Z_COP.l;
ranges('OF1max')=Z.l;

E=ranges('OF1max');
loop(counter,
E=(ranges('OF2max')-ranges('OF2min'))*(ord(counter)-1)/(card(counter)-1)+ranges('OF2min');
*E=ranges('OF2max')-(ranges('OF2max')-ranges('OF2min'))*(ord(counter)-1)/(card(counter)-1);
Z_CO2.up=E;
*Z_COP.up=E;
option resLim = 4800
option optcr = 0.003
Solve DHC_with_shared_inv using mip minimizing Z ;

SCOP(p_nDC)$(P_tot.l(p_nDC)>0)= Q_tot.l(p_nDC)/(P_tot.l(p_nDC));
SCOP('DC')$(P_tot.l('DC')>0)= Q_tot.l('DC')/(P_tot.l('DC')*2);
SCOP(p)$(P_tot.l(p)=0)= 0;

SCOP_tot$(SCOP_coef=0) = 0;
SCOP_tot$(SCOP_coef>0) = sum(t,d(t))/conv/SCOP_coef;

* Cost of HP production
C_prod_HP_p(p)$(Q_tot.l(p)>0) = ((C_inv_tot.l(p)+C_DH.l(p))*a_c_inv  + C_m.l(p) + C_el_p.l(p) )/(Q_tot.l(p)*conv)+C_inv_stor_p(p)+C_prod_DH_net_p(p);
C_prod_HP_p(p)$(Q_tot.l(p)=0) = 0;
C_prod_HP = ((C_inv_tot_all.l+C_inv_stor.l+C_inv_DH_net+C_DH_pipe.l)*a_c_inv + C_m_tot.l + C_el_tot.l)/sum(t,d(t));
CO2_HP_p(p)$(Q_tot.l(p)>0) = sum(t,(P_hs.l(p,t))*CO2(t)/conv)/Q_tot.l(p) ;
CO2_HP_p(p)$(Q_tot.l(p)=0)= 0;
*CO2_HP_p(p) = sum(t,(P_hs.l(p,t)+P_aux_el_t.l(p,t))*CO2(t)/conv) ;
CO2_HP = sum((p,t),(P_hs.l(p,t))*CO2(t))/sum(t,d(t));

SCOP_c(p_nDC)$(P_tot_c.l(p_nDC)>0)= (Q_tot_c.l(p_nDC)+Q_free_tot.l(p_nDC))/(P_tot_c.l(p_nDC)+P_free_tot.l(p_nDC) );
SCOP_c('DC')$(P_tot_c.l('DC')>0)= (Q_tot_c.l('DC')+Q_free_tot.l('DC'))/(P_tot_c.l('DC')*2+P_free_tot.l('DC') );
SCOP_c(p)$(P_tot_c.l(p)=0)= 0;

SCOP_tot_c$(SCOP_coef1=0) = 0;
SCOP_tot_c$(SCOP_coef1>0) = sum(t,d_DC(t))/conv/SCOP_coef1;
SCOP_sys$(SCOP_coef2=0) = 0;
SCOP_sys$(SCOP_coef2>0)=  sum(t,d_DC(t)+d(t))/conv/SCOP_coef2;

* Cost of HP production
C_prod_HP_p_c(p)$(Q_tot_c.l(p)>0) = ((C_inv_tot_c.l(p)+C_DH_c.l(p)+C_inv_free.l(p))*a_c_inv  + C_m_c.l(p) + C_el_p_c.l(p) )/((Q_tot_c.l(p)+Q_free_tot.l(p))*conv)+C_inv_stor_p_c(p);
C_prod_HP_p_c(p)$(Q_tot_c.l(p)=0) = 0;
C_prod_HP_c = ((C_inv_tot_all_c.l+C_inv_stor_c.l+C_DH_pipe_c.l+C_inv_free_tot.l)*a_c_inv + C_m_tot_c.l + C_el_tot_c.l )/sum(t,d_DC(t));
CO2_HP_p_c(p)$(Q_tot_c.l(p)>0) = sum(t,(P_hs_c.l(p,t)+P_free.l(p,t))*CO2(t)/conv)/(Q_tot_c.l(p)+Q_free_tot.l(p)) ;
CO2_HP_p_c(p)$(Q_tot_c.l(p)=0)= 0;
CO2_HP_c = sum((p,t),(P_hs_c.l(p,t)+P_free.l(p,t))*CO2(t))/sum(t,d_DC(t));

report(counter,'Z')=Z.l;
report(counter,'Z_CO2')=Z_CO2.l;
report(counter,'Z_COP')=Z_COP.l;
report(counter,'E')=E;
report(counter,'St_V')=St_V.l;
report(counter,'St_c')=St_c.l;
report(counter,'St_V_c')=St_V_c.l;
report(counter,'St_c_c')=St_c_c.l;
report(counter,'LCOH')=C_prod_HP;
report(counter,'LCOC')=C_prod_HP_c;
report(counter,'CO2_H')=CO2_HP;
report(counter,'CO2_C')=CO2_HP_c;
report(counter,'SCOP')=SCOP_tot;
report(counter,'SCOP_c')=SCOP_tot_c;
report(counter,'SCOP_sys')=SCOP_sys;

report_p(counter,'Q_sink_d',p)=Q_sink_d.l(p);
report_p(counter,'Q_sink_d_c',p)=Q_sink_d_c.l(p);
report_p(counter,'Q_free_d',p)=Q_free_d.l(p);
report_p(counter,'LCOH_p',p)=C_prod_HP_p(p);
report_p(counter,'LCOC_p',p)=C_prod_HP_p_c(p);
report_p(counter,'CO2_H_p',p)=CO2_HP_p(p);
report_p(counter,'CO2_C_p',p)=CO2_HP_p_c(p);
report_p(counter,'SCOP_p',p)=SCOP(p);
report_p(counter,'SCOP_p_c',p)=SCOP_c(p);

);
Display report,report_p;

$offtext

***************************************************************
*$offtext

*EXECUTE_UNLOAD 'GAMS_output_DHC_pareto';


