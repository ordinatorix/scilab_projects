// Project 2: Ramjet Design

// Close everything-------------------------------------------------------
clear; clc;close;close;close;close;close; close;

// Independent Parameters-------------------------------------------------
T_a         = 216.65; // (K)
P_a         = 18822.653374; //19193.37995; // (Pa)
rho_a       = 0.30266; // (kg/m^3)
r_rho       = 0.24708;
M_a         = 4; //[1.25:0.25:4];
Thrust      = 44482.22; // (N)
gamma_a     = 1.4;
gamma_c     = 1.4; //1.33;//   ASSUMING GAMMA IS CONSTANT
gamma_e     = 1.4; //1.36;
R           = 287; // (J/kgK)
FA_st       = 0.06; 
To_4        = 2300; // (K)
P_e         = P_a;
cp_a        = (gamma_a/(gamma_a-1))*R;
cp_c        = (gamma_c/(gamma_c-1))*R;
cp_e        = (gamma_e/(gamma_e-1))*R;
M_e_min     = M_a;
M_5         = 1;

// Variable parameters
A_e         = [0.01:0.01:1];
M_e         = [4:0.00201:4.2];//[4:0.0135:5.34]; //[4:0.0201:6]; // SINCE GAMMA IS CONSTANT, M_e = M_a
r_b         = [0.81:0.01:0.99];
r_n         = 1;

// Dependent Parameters---------------------------------------------------
Po_6          = zeros(1:100);
A_star        = zeros(100,21);
A_5           = zeros(100,21);
A_1           = zeros(100,21);
A_2           = zeros(100,21);
A_4           = zeros(100,21);
m_dot_a       = zeros(100,21);
m_dot_e       = zeros(100,21);
m_dot_f       = zeros(100,21);

// Ambient Air Conditions-------------------------------------------------
Po_a        = P_a*(1+((gamma_a-1)/2.0)*M_a.^2).^((gamma_a)/(gamma_a-1));
To_a        = T_a.*(1+(gamma_a-1)./2.0.*M_a.^2);
u_a         = M_a.*sqrt(gamma_a*R*T_a);

// Inlet (a-1)   ---------------------------------------------------------
r_d         = 1-0.075.*(M_a-1).^1.35; // Based on complex shock structure
Po_2        = Po_a.*r_d;;
To_2        = To_a;

// Nozzle Exhaust Conditions----------------------------------------------

for i = 1:100,
    for j = 1:100,
        A_star(i,j)      = A_e(i)./((1.0./M_e(j)).*((2.0./(gamma_e+1)).*(1+((gamma_e-1)./2.0).*M_e(j).^2)).^((gamma_e+1)./(2.0.*(gamma_e-1))));
        A_5         = A_star;
        A_1         = A_5;
        A_2         = A_1;
    end
end

Po_6        = P_e.*(1+((gamma_e-1)./2).*M_e.^2).^(gamma_e./(gamma_e-1));
P_6         = P_e;
To_6        = To_4;
p_r_6       = Po_6./P_6;
T_6         = To_6./(p_r_6.^((gamma_e-1)/gamma_e));
u_e         = M_e.*sqrt(gamma_e.*R.*T_6);

//Calculate Mass Flow Rate
m_dot_a     = u_a.*(P_a/(R*T_a)).*A_1;
m_dot_e     = (P_6./(R.*T_6)).*u_e.*A_e;

for i = 1:100,
    for j = 1:100,
        m_dot_ei(i,j) = m_dot_e(j);
    end
end

for i = 1:100,
    for j = 1:100,
        M_ei(i,j) = M_e(j);
    end
end

// Throat Conditions------------------------------------------------------
M_5         = 1;
Po_5        = Po_6;
To_5        = To_4;
P_5         = Po_5./(1+((gamma_e-1)/2).*M_5.^2).^(gamma_e./(gamma_e-1));
T_5         = To_5./(1+((gamma_e-1)/2).*M_5.^2);
u_5         = M_5.*sqrt(gamma_e.*R.*T_5);

// Burner Exit Conditions-------------------------------------------------
M_4         = [0.1:0.00891:0.99];
Po_4        = Po_6; 
m_dot_f     = m_dot_ei-m_dot_a;
FA          = (m_dot_ei./m_dot_a)-1;
for i = 1:100,
    for j = 1:100,
        if FA(i,j) < 0 then
            FA(i,j)= 0.0000001;
        elseif FA(i,j) > 0.06 then
            FA(i,j) = 0.06;
        end
    end
end

eq_ratio    = FA./FA_st;
for i = 1:100,
    for j = 1:100,
        if eq_ratio(i,j) > 1  then
            messagebox("Equivalence Ratio Impossible!","!ERROR!","error");
            disp(eq_ratio);
            abort;
        end
    end
end


for i = 1:100,
    for j = 1:100,
        u_ei(i,j) = u_e(j);
    end
end

Thrust_check = m_dot_ei.*u_ei-m_dot_a.*u_a;

for i = 1:100,
    for j = 1:100,
        if Thrust_check(i,j) < 0 then
            Thrust_check(i,j)=0;
        elseif Thrust_check(i,j) > Thrust  then
            T =1;
        end
    end
end

if  T == 1 then
    messagebox("Maximum Thrust Exceeded!","!WARNING!","warning");
end
for i = 1:100,
    for j = 1:100,
        M_4i(i,j) = M_4(j);
    end
end

A_4         = A_5.*((1.0./M_4i).*((2.0./(gamma_c+1)).*(1+((gamma_c-1)./2.0).*M_4i.^2)).^((gamma_c+1)./(2.0.*(gamma_c-1))));
A_3         = A_4;

if A_3 > 10.*A_1 then
    messagebox("Maximum Combustor Area Exceeded!","!WARNING!","warning");
end

TSFC           = FA./((1+FA).*u_ei-u_a);
h_c         = cp_c.*To_2.*((((To_4./To_2)-1)./FA)+(To_4./To_2));

f0 = scf(0); //creates figure with id==0 and make it the current one
f1 = scf(1); //creates figure with id==1 and make it the current one
f2 = scf(2);

// plot Thrust VS M_e VS TSFC
scf(f0);
subplot(111);
param3d(Thrust_check.*10^(-5),TSFC.*10^5,M_ei);
xtitle("TSFC VS M_e VS Thrust","Thrust","TSFC","M_e");

scf(f1);
subplot(111);
param3d(A_5,Thrust_check.*10^(-6),TSFC.*10^4);
xtitle("TSFC VS Area 5 VS Thrust","A_5","Thrust","TSFC");

scf(f2);
subplot(111);
plot(A_5,TSFC);
xtitle("TSFC VS Area 5","A_5","TSFC");




//col 83
    //M_e = 4.607
    //A_e = 0.84
    //A_5 = 0.08

