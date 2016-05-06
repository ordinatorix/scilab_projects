// Hw 8-------------------------------------------------------------------

// Close and clear all
clc; clear; close; close; close; close; close; close; close; close; close;
close;

// Independant Parameters------------------------------------------------
M_a         = 0.8;
P_a         = 26500; // (Pa)
T_a         = 223.252; // (K)
h_c         = 43000000; // (j/kg)
gamma_a     = 1.4;
gamma_c     = 1.35;
gamma_e     = gamma_c;
R           = 287;
r_f         = 2;
r_c         = 24;
r_b         = 0.97;
FA_st       = 0.06;
To_4_max    = 1500; // (K)
eff_diff    = 0.94;
eff_comp    = 0.87;
eff_fan     = 0.92;
eff_burn    = 0.98;
eff_turb    = 0.85;
eff_c_nzl   = 0.97;
eff_f_nzl   = 0.98;
cp_a        = (gamma_a/(gamma_a-1))*R;
cp_c        = (gamma_c/(gamma_c-1))*R;
cp_         = (gamma_e/(gamma_e-1))*R;
P_e         = P_a;

// Design Variable--------------------------------------------------------
B           = [0.1:0.1:7.2];
//B = [0.1:0.1:20];

// Dependant Paramenters--------------------------------------------------
W_t_out     = zeros(1,72);
To_5_r      = zeros(1,72);
To_5_i      = zeros(1,72);
Po_5        = zeros(1,72);
To_7_r      = zeros(1,72);
Po_7        = zeros(1,72);
P_7         = zeros(1,72);
T_7_i       = zeros(1,72);
T_7_r       = zeros(1,72);
u_c_e       = zeros(1,72);
M_c_e       = zeros(1,72);
Thrust      = zeros(1,72);
I           = zeros(1,72);
TSCF        = zeros(1,72);
eff_propul  = zeros(1,72);
eff_therm   = zeros(1,72);
ell_overall = zeros(1,72);


// Ambient Conditions-----------------------------------------------------
u_a         = M_a*sqrt(gamma_a*R*T_a);
Po_a        = P_a*(1+(gamma_a-1)/2*M_a^2)^(gamma_a/(gamma_a-1));
To_a        = T_a*(1+(gamma_a-1)/2*M_a^2);

// Inlet/Diffuser Conditions----------------------------------------------
To_2_r      = To_a;
To_2_i      = eff_diff*(To_2_r-T_a)+T_a;
Po_2        = P_a*(To_2_i/T_a)^(gamma_a/(gamma_a-1));

// Fan Outlet Conditions---------------------------------------------------------
Po_8        = r_f*Po_2;
To_8_i      = To_2_r*(r_f)^((gamma_a-1)/gamma_a);
To_8_r      = (To_8_i-To_2_r)/eff_fan+To_2_r;
W_f_in      = B.*cp_a.*(To_8_r-To_2_r);

// Compressor Conditions--------------------------------------------------
Po_3        = r_c*Po_2;
To_3_i      = To_2_r*r_c^((gamma_a-1)/gamma_a);
To_3_r      = (To_3_i-To_2_r)/eff_comp+To_2_r;
W_c_in      = cp_a*(To_3_r-To_2_r);

// Burner Conditions------------------------------------------------------
FA          = ((To_4_max/To_3_r-1)/((eff_burn*h_c)/(cp_c*To_3_r)-(To_4_max/To_3_r)));
eq_ratio    = FA/FA_st;
To_4_r      = To_4_max;
Po_4        = r_b*Po_3;

// Turbine Conditions-----------------------------------------------------
W_t_out     = W_c_in+W_f_in;
To_5_r      = To_4_r-W_t_out./((1+FA).*cp_c);
To_5_i      = To_4_r-(To_4_r-To_5_r)./eff_turb;
Po_5        = Po_4.*(To_5_i./To_4_r).^(gamma_c./(gamma_c-1));

// Core Nozzle Conditions------------------------------------------------------
To_7_r      = To_5_r;
Po_7        = Po_5;
P_7         = P_e;
T_7_i       = To_7_r.*(P_7./Po_7).^((gamma_e-1)./gamma_e);
T_7_r       = To_5_r-((To_5_r-T_7_i).*eff_c_nzl);
u_c_e       = sqrt(2.0.*eff_c_nzl.*(gamma_e./(gamma_e-1)).*R.*To_5_r.*(1-(P_e./Po_5).^((gamma_e-1)./(gamma_e))));
M_c_e       = u_c_e./sqrt(gamma_e.*R.*T_7_r);
A_ratio_c   = 1.0./M_c_e.*((2.0./(gamma_e+1)).*(1+((gamma_e-1)./2.0).*M_c_e.^2)).^((gamma_e+1)./(2.0.*(gamma_e-1)));

// Fan Nozzle Conditions--------------------------------------------------
To_9_r      = To_8_r;
T_9_i       = To_8_r./((Po_8./P_e).^((gamma_a-1)./gamma_a));
T_9_r       = To_8_r-(eff_f_nzl.*(To_8_r-T_9_i));
u_f_e       = sqrt(2.0.*eff_f_nzl.*(gamma_a./(gamma_a-1)).*R.*To_8_r.*(1-(P_a./Po_8).^((gamma_a-1)./(gamma_a))));
M_f_e       = u_f_e./sqrt(gamma_a.*R.*T_9_r);
A_ratio_f   = 1.0./M_f_e.*((2.0./(gamma_e+1)).*(1+((gamma_e-1)./2.0).*M_f_e.^2)).^((gamma_e+1)./(2.0.*(gamma_e-1)));

// Engine Performance-----------------------------------------------------
Thrust      = B.*(u_f_e-u_a)+((1+FA).*u_c_e-u_a);
I           = (1+FA).*u_c_e+B.*u_f_e-(1+B).*u_a; 
TSFC        = FA./I;
W_dot_P     = Thrust.*u_a;
m_dot       = 1;
P_avail     = m_dot.*((1+FA).*((u_c_e.^2)./2)+B.*((u_f_e.^2)./2)-(B+1).*(u_a.^2)./2);
P_in        = FA.*m_dot.*h_c;
eff_propul  =(I.*u_a)./((1+FA).*((u_c_e.^2)./2)+B.*((u_f_e.^2)./2)-(B+1).*(u_a.^2)./2);
eff_therm   = P_avail./P_in;
eff_overall = eff_therm.*eff_propul;

// Plot Graphs to be Compared---------------------------------------------
f0 = scf(5); //creates figure with id==0 and make it the current one
f1 = scf(6); //creates figure with id==1 and make it the current one
f2 = scf(7);
f3 = scf(8);
f4 = scf(9);

// Plot Specific Thrust---------------------------------------------------
scf(f0);
subplot(311);
plot(B,I, 'r');
xtitle("Specific Thurst VS Bipass Ratio","B","Specific Thurst(N*s/Kg)");

// plot TSFC--------------------------------------------------------------
subplot(312);
plot(B,TSFC, 'g');
xtitle("TSFC VS Bipass Ratio","B","TSFC(Kg/N*s)");

// plot area ratio--------------------------------------------------------
subplot(313);
plot(B,A_ratio_c, 'r');
plot(B,A_ratio_f, 'b');
xtitle("Area Ratio VS Bipass Ratio","B","Area Ratio");


// Plot Temperature Across TurboJet---------------------------------------
scf(f1);
subplot(111);
plot(B,To_5_r, 'c');
plot(B,T_7_r, 'm');
xtitle("Stagnation Temperature VS Bipass Ratio","B","Temperature (K)");
leg=legend(['To_5_r';'T_7_r'],[3]);

// Plot Temperature Across TurboJet---------------------------------------
scf(f2);
subplot(111);
plot(B,Po_5, 'm');
xtitle("Stagnation Pressure VS Bipass Ratio","B","Pressure (Pa)");
leg=legend(['Po_5'],[4]);

// plot efficiencies------------------------------------------------------
scf(f3); 
subplot(111);
plot(B,eff_therm, 'r');
plot(B,eff_propul, 'b');
plot(B,eff_overall, 'g');
xtitle("Efficiencies VS Bipass Ratio","B","Effeciency");
leg=legend(['Thermal Efficiency';'Propulsive Efficiency';'Overall Efficiency'],[2]);

// Plot Exit Mach number--------------------------------------------------
scf(f4);
subplot(111);
plot(B,M_c_e, 'r');
plot(B,M_f_e, 'r');
xtitle("Core & Fan Exit Mach Number VS Bipass Ratio","B","M");

// optimal bypass ratio at 5.9










