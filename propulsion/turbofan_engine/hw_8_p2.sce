// Hw 8-------------------------------------------------------------------

// Close and clear all
clc; clear; close; close; close; close; close; close; close; close; close;
close;

// Using the optimal bypass ratio of 5.9 found in problem 1---------------

// Independant Parameters-------------------------------------------------
M_a         = 0.8;
P_a         = 26500; // (Pa)
T_a         = 223.252; // (K)
h_c         = 43000000; // (j/kg)
gamma_a     = 1.4;
gamma_c     = 1.35;
gamma_e     = gamma_c;
R           = 287;
r_f         = [1.5:0.1:2.2];
r_c         = [20:0.1:28]; 
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
B           = 5.9;

// Dependant Paramenters--------------------------------------------------

Po_8        = zeros(1,8);
To_8_r      = zeros(1,8);
To_8_i      = zeros(1,8);
W_f_in      = zeros(1,8);


FA          = zeros(1,81);
eq_ratio    = zeros(1,81);

W_t_out     = zeros(8,81);
To_5_r      = zeros(8,81);
To_5_i      = zeros(8,81);
Po_5        = zeros(8,81);
To_7_r      = zeros(8,81);
Po_7        = zeros(8,81);
P_7         = zeros(8,81);
T_7_i       = zeros(8,81);
T_7_r       = zeros(8,81);
u_c_e       = zeros(8,81);
M_c_e       = zeros(8,81);
A_ratio_c   = zeros(8,81);
Thrust      = zeros(8,81);
I           = zeros(8,81);
TSCF        = zeros(8,81);
P_in        = zeros(8,81);
P_avail     = zeros(8,81);
W_dot_P     = zeros(8,81);
eff_propul  = zeros(8,81);
eff_therm   = zeros(8,81);
ell_overall = zeros(8,81);


// Ambient Conditions-----------------------------------------------------
u_a         = M_a*sqrt(gamma_a*R*T_a);
Po_a        = P_a*(1+(gamma_a-1)/2*M_a^2)^(gamma_a/(gamma_a-1));
To_a        = T_a*(1+(gamma_a-1)/2*M_a^2);

// Inlet/Diffuser Conditions----------------------------------------------
To_2_r      = To_a;
To_2_i      = eff_diff*(To_2_r-T_a)+T_a;
Po_2        = P_a*(To_2_i/T_a)^(gamma_a/(gamma_a-1));

// Fan Outlet Conditions---------------------------------------------------------
Po_8        = r_f.*Po_2;
To_8_i      = To_2_r.*(r_f).^((gamma_a-1)./gamma_a);
To_8_r      = (To_8_i-To_2_r)./eff_fan+To_2_r;
W_f_in      = B.*cp_a.*(To_8_r-To_2_r);

// Compressor Conditions--------------------------------------------------
Po_3        = r_c.*Po_2;
To_3_i      = To_2_r.*r_c.^((gamma_a-1)./gamma_a);
To_3_r      = (To_3_i-To_2_r)./eff_comp+To_2_r;
W_c_in      = cp_a.*(To_3_r-To_2_r);

// Burner Conditions------------------------------------------------------
FA          = ((To_4_max./To_3_r-1)./((eff_burn.*h_c)./(cp_c.*To_3_r)-(To_4_max./To_3_r)));
eq_ratio    = FA./FA_st;
To_4_r      = To_4_max;
Po_4        = r_b.*Po_3;

// Turbine Conditions-----------------------------------------------------
for i = 1:8,
    for j = 1:81,
        W_t_out(i,j)= W_c_in(j)+W_f_in(i);
        To_5_r(i,j) = To_4_r-W_t_out(i,j)./((1+FA(j)).*cp_c);
        To_5_i(i,j) = To_4_r-(To_4_r-To_5_r(i,j))./eff_turb;
        Po_5(i,j)   = Po_4(j).*(To_5_i(i,j)./To_4_r).^(gamma_c./(gamma_c-1));
    end
end

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
for i = 1:8,
    for j = 1:81,
        Thrust(i,j)      = B.*(u_f_e(i)-u_a)+((1+FA(j)).*u_c_e(i,j)-u_a);
        I(i,j)           = (1+FA(j)).*u_c_e(i,j)+B.*u_f_e(i)-(1+B).*u_a; 
        TSFC(i,j)        = FA(j)./I(i,j);
        W_dot_P(i,j)     = Thrust(i,j).*u_a;
        m_dot            = 1;
        P_avail(i,j)     = m_dot.*((1+FA(j)).*((u_c_e(i,j).^2)./2)+B.*((u_f_e(i).^2)./2)-(B+1).*(u_a.^2)./2);
        P_in(i,j)        = FA(j).*m_dot.*h_c;
        eff_propul(i,j)  =(I(i,j).*u_a)./((1+FA(j)).*((u_c_e(i,j).^2)./2)+B.*((u_f_e(i).^2)./2)-(B+1).*(u_a.^2)./2);
        eff_therm(i,j)   = P_avail(i,j)./P_in(i,j);
        eff_overall(i,j) = eff_therm(i,j).*eff_propul(i,j);
    end
end

// Plot Graphs to be Compared---------------------------------------------
f0 = scf(0); //creates figure with id==0 and make it the current one
f1 = scf(1); //creates figure with id==1 and make it the current one
f2 = scf(2);
f3 = scf(3);
f4 = scf(4);
f5 = scf(5);
f6 = scf(6);
f7 = scf(7);
f8 = scf(8);

// Plot Specific Thrust---------------------------------------------------
scf(f0);
subplot(311);
plot3d(r_f,r_c,I, 'r');
xtitle("Specific Thurst VS Fan & Compressor Pressure Ratio at 5.9 Bypass ratio","r_f","r_c","Specific Thurst(N*s/Kg)");

// plot TSFC--------------------------------------------------------------
subplot(312);
plot3d(r_f,r_c,TSFC, 'g');
xtitle("TSFC VS Fan & Compressor Pressure Ratio at 5.9 Bypass ratio","r_f","r_c","TSFC(Kg/N*s)");

// plot area ratio--------------------------------------------------------
subplot(313);
plot3d(r_f,r_c,A_ratio_c, 'r');
xtitle("Area Ratio VS Fan & Compressor Pressure Ratio at 5.9 Bypass ratio","r_f","r_c","Area Ratio");

scf(f1)
subplot(111);
plot(A_ratio_f,r_f, 'b');
xtitle("Area Ratio VS Fan & Compressor Pressure Ratio at 5.9 Bypass ratio","Area Ratio","r_f");

// Plot Temperature Across TurboJet---------------------------------------
scf(f2);
subplot(111);
plot3d(r_f,r_c,To_5_r, 'c');
plot3d(r_f,r_c,T_7_r, 'm');
xtitle("Stagnation Temperature VS Fan & Compressor Pressure Ratio at 5.9 Bypass ratio","r_f","r_c","Temperature (K)");


// Plot Temperature Across TurboJet---------------------------------------
scf(f3);
subplot(111);
plot3d(r_f,r_c,Po_5, 'm');
xtitle("Turbine Exhaust Stagnation Pressure VS Fan & Compressor Pressure Ratio at 5.9 Bypass ratio","r_f","r_c","Pressure (Pa)");


// plot efficiencies------------------------------------------------------
scf(f4); 
subplot(111);
plot3d(r_f,r_c,eff_overall.*100);
xtitle("Overall Efficiencies VS Fan & Compressor Pressure Ratio at 5.9 Bypass ratio","r_f","r_c","Effeciency(%)");


scf(f5); 
subplot(111);
plot3d(r_f,r_c,eff_therm.*100, 'r');
xtitle("Thermal Efficiencies VS Fan & Compressor Pressure Ratio at 5.9 Bypass ratio","r_f","r_c","Effeciency(%)");

scf(f6); 
subplot(111);
plot3d(r_f,r_c,eff_propul.*100, 'b');
xtitle("Propulsive Efficiencies VS Fan & Compressor Pressure Ratio at 5.9 Bypass ratio","r_f","r_c","Effeciency(%)");

// Plot Exit Mach number--------------------------------------------------
scf(f7);
subplot(111);
plot3d(r_f,r_c,M_c_e);
xtitle("Core Exit Mach Number VS Fan & Compressor Pressure Ratio at 5.9 Bypass ratio","r_f","r_c","M");

scf(f8);
subplot(111);
plot(M_f_e,r_f, 'r');
xtitle("Fan Exit Mach Number VS Fan Pressure Ratio at 5.9 Bypass ratio","M","r_f");

// For the folowing values, the overall efficiency of the turbofan reachs a maximum
// r_f = 1.9
// r_c = 28










