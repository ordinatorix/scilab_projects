// clear previous calculations
close; clc; clear; close;

// Independent Parameters
gamma_a = 1.4;
gamma_e = 1.33;
R = 287;
T_a = 223.252;
P_a = 16500;
h_c = 43000000;
FA_st = 0.06;
T_max = 2600;
cp_a = (gamma_a/(gamma_a-1)*R);
cp_e = (gamma_e/(gamma_e-1)*R);
M = [1:0.25:6];

//setting up empty matrices
To_4           = zeros(1,21);
T_e            = zeros(1,21);
u_e            = zeros(1,21);
I              = zeros(1,21);
TSFC           = zeros(1,21);
thermal_eff    = zeros(1,21);
propulsive_eff = zeros(1,21);
overall_eff    = zeros(1,21);
M_e            = zeros(1,21);

// Dependent Parameters---------------------------------------

// Find flow velocity
u= M*sqrt(gamma_a*R*T_a);

// Inlet/Diffuser
To_2 = T_a*(1+(gamma_a-1)/2 *M.^2);

// Calculate fuel to air ratio
FA = ((T_max./To_2)-1)./((h_c./(cp_e.*To_2))-(T_max./To_2));

// Find equivalence ratio and make sure it is less than or equal to unity.
eq_ratio = FA./FA_st;
for i = 1:21,
    if eq_ratio(i) <= 1 then
        To_4(i) = T_max;
    else
    // Calculate To_4 at FA_st
        To_4(i) = ((FA_st.*(h_c./(cp_e.*To_2(i)))+1)./(FA_st+1)).*To_2(i);
    end,
end

// Calculate T_e
T_e = (T_a./To_2).*To_4;

// From T_e calculate u_e
u_e = sqrt(2*cp_e.*(To_4-T_e));

// Calculate exit Mach number
M_e = u_e./sqrt(gamma_e.*R.*T_e);

// Calculate Specific thrust; TSFC; Thermal Efficiency; Propulsive Efficiency &

// Overall Efficiency
for  i = 1:21,
    if eq_ratio(i) > 1  then
        FA(i) = FA_st;
    end
    I(i)              = ((1+FA(i)).*u_e(i)-u(i));
    TSFC(i)           = FA(i)./((1+FA(i)).*u_e(i)-u(i));
    thermal_eff(i)    = ((1+FA(i)).*u_e(i).^2-u(i).^2)./(2.0.*FA(i).*h_c);
// I get a thermal efficiency of more than 1. see with TA for correctness. If 
//cannot be greater than 1, implement if statement.
//    if thermal_eff(i) >1 then
//        thermal_eff(i) = 1;
//    end
    propulsive_eff(i) = (I(i).*u(i))./((1+FA(i)).*(u_e(i).^2.0./2)-(u(i).^2.0./2));
    overall_eff(i)    = thermal_eff(i)*propulsive_eff(i);
end

//Calculate the area ratio
A_ratio = 1.0./M_e.*((2/(gamma_e+1)).*(1+(gamma_e-1)./2.0.*M.^2)).^((gamma_e+1)./(2.0.*(gamma_e-1)));




//Plots-------------------------------------------------
f0 = scf(0); //creates figure with id==4 and make it the current one
f1 = scf(1); //creates figure with id==0 and make it the current one
//f2 = scf(2);


// plot specific thrust
scf(f0);
subplot(411);
plot(M,I, 'r');
xtitle("Specific Thurst VS Mach Number","Mach Number","Specific Thurst(N*s/Kg)");
a=get("current_axes")
a.data_bounds=[1,500;6,1220];

// plot TSFC
subplot(412);
plot(M,TSFC, 'g');
xtitle("TFSC VS Mach Number","Mach Number","TFSC(Kg/N*s)");

// plot combustor exit temperature
subplot(413);
plot(M,To_4);
xtitle("To_4 VS Mach Number","Mach Number","To_4 (K)");
a=get("current_axes")
a.data_bounds=[1,2300;6,2800];

// plot area ratio
subplot(414);
plot(M,A_ratio, 'r');
xtitle("Area Ratio VS Mach Number","Mach Number","Area Ratio");

// plot efficiencies
scf(f1); // set first created figure as current one
subplot(111);
plot(M,thermal_eff, 'r');
plot(M,propulsive_eff);
plot(M,overall_eff, 'g');
xtitle("Efficiencies VS Mach Number","Mach Number","Effeciency");
leg=legend(['Thermal Efficiency';'Propulsive Efficiency';'Overall Efficiency'],[4]);
