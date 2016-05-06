// clear previous calculations
close; clc; clear; close; close; close

// Independent Parameters
gamma_a = 1.4;
R = 287;
T_a = 220;
P_a = 8500;
h_c = 43000000;
To_4 = 2540;
cp_a = (gamma_a/(gamma_a-1)*R);
M = [1:0.25:6];
combustion_eff = [0.01:0.01:1];
r_n = [0.01:0.01:1];

//setting up empty matrices
To_2           = zeros(1,100);
T_e            = zeros(1,100);
u_e            = zeros(1,100);
I              = zeros(1,100);
TSFC           = zeros(1,100);
I_comb_eff     = zeros(21,100);
I_r_n_eff      = zeros(21,100);
TSFC_comb_eff  = zeros(21,100);
TSFC_r_n_eff   = zeros(21,100);


// Dependent Parameters---------------------------------------

// Find Flow Velocity
u    = M*sqrt(gamma_a*R*T_a);

// Find Diffuser/Combustor Stagnation Temperature
To_2 = T_a.*(1+(gamma_a-1)./2.0.*M.^2);

// Find Nozzle Exit Temperature
T_e  = T_a./To_2.*To_4;

// Find Nozzle Exit Velocity
u_e  = sqrt(2.0.*cp_a*(To_4-T_e));

// Calculate Specific Thrust
I    = u_e-u;

// Calculate TSFC
TSFC = 1.0./I;

//Calculate d(I)/comb_eff; d(I)/r_n; ---
for i = 1:21,
    for j = 1:100,
        I_comb_eff(i,j)    = I(i)./combustion_eff(j);
        I_r_n_eff(i,j)     = I(i)./r_n(j);
        TSFC_comb_eff(i,j) = TSFC(i)./combustion_eff(j);
        TSFC_r_n_eff(i,j)  = TSFC(i)./r_n(j);
    end,
end


//plot------------------------------------------------------------------------
f0 = scf(0); //creates figure with id==0 and make it the current one
f1 = scf(1); //creates figure with id==1 and make it the current one
f2 = scf(2);
f3 = scf(3);

// plot d(I)/Comb_eff VS Mach Number
scf(f0);
subplot(111);
plot(M,I_comb_eff, 'g');
xtitle("d(I)/Comb_eff VS Mach Number","Mach Number","d(I)/Comb_eff");
//a=get("current_axes")
//a.data_bounds=[1,500;6,1220];

// plot d(I)/Nozzle pressure ratio VS Mach Number
scf(f1);
subplot(111);
plot(M,I_r_n_eff, 'r');
xtitle("d(I)/r_n VS Mach Number","Mach Number","d(I)/r_n");

// plot d(TSFC)/Comb_eff VS Mach Number
scf(f2);
subplot(111);
plot(M,TSFC_comb_eff, 'g');
xtitle("d(TSFC)/Comb_eff VS Mach Number","Mach Number","d(TSFC)/Comb_eff");

// plot d(TSFC)/Nozzle Pressure Ratio VS Mach Number
scf(f3);
subplot(111);
plot(M,TSFC_r_n_eff, 'r');
xtitle("d(TSFC)/r_n VS Mach Number","Mach Number","d(TSFC)/r_n");

