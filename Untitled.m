%%%%%%%  HW4 Q3 %%%%%%
clc; clear;
%%% given data %%%

%% all units in SI%%
Gna = 1.2;                  %[mS/mm^2]
Gk = 0.36;                  %[mS/mm^2]
Gl = 0.003;                 %[mS/mm^2]
Ena = 50;                   %[mV]
Ek = -77;                   %[mV]
El = -54.4;                 %[mV]
C = 0.01;                   %[uF/mm^2]
V1_0 = -60;                 %[mV]
V2_0 = -60;                 %[mV]
V3_0 = -60;                 %[mV]
V4_0 = -60;                 %[mV]
dt = 0.01;                  %[msec]
L = 10;                     %[um]
radius1 = 1;                %[um]
radius2 = 1;                %[um]
radius3 = 2;                %[um]
radius4 = 0.2;              %[um]
r_L = 1;                    %[kohm*mm] 
r_m = 1;                    %[Mohm*mm^2]

%% set wanted vectors & initial parameters %%
V1 = zeros(1,100000);
V2 = zeros(1,100000);
V3 = zeros(1,100000);
V4 = zeros(1,100000);

n1 = zeros(1,100000);
m1 = zeros(1,100000);
h1 = zeros(1,100000);

n2 = zeros(1,100000);
m2 = zeros(1,100000);
h2 = zeros(1,100000);

n3 = zeros(1,100000);
m3 = zeros(1,100000);
h3 = zeros(1,100000);

n4 = zeros(1,100000);
m4 = zeros(1,100000);
h4 = zeros(1,100000);

alpha_n1 = zeros(1,100000);
alpha_m1 = zeros(1,100000);
alpha_h1 = zeros(1,100000);
beta_n1 = zeros(1,100000);
beta_m1 = zeros(1,100000);
beta_h1 = zeros(1,100000);

alpha_n2 = zeros(1,100000);
alpha_m2 = zeros(1,100000);
alpha_h2 = zeros(1,100000);
beta_n2 = zeros(1,100000);
beta_m2 = zeros(1,100000);
beta_h2 = zeros(1,100000);

alpha_n3 = zeros(1,100000);
alpha_m3 = zeros(1,100000);
alpha_h3 = zeros(1,100000);
beta_n3 = zeros(1,100000);
beta_m3 = zeros(1,100000);
beta_3 = zeros(1,100000);

alpha_n4 = zeros(1,100000);
alpha_m4 = zeros(1,100000);
alpha_h4 = zeros(1,100000);
beta_n4 = zeros(1,100000);
beta_m4 = zeros(1,100000);
beta_h4 = zeros(1,100000);

tau_n1 = zeros(1,100000);
tau_m1 = zeros(1,100000);
tau_h1 = zeros(1,100000);
n_eq1 = zeros(1,100000);
m_eq1 = zeros(1,100000);
h_eq1 = zeros(1,100000);

tau_n2 = zeros(1,100000);
tau_m2 = zeros(1,100000);
tau_h2 = zeros(1,100000);
n_eq2 = zeros(1,100000);
m_eq2 = zeros(1,100000);
h_eq2 = zeros(1,100000);

tau_n3 = zeros(1,100000);
tau_m3 = zeros(1,100000);
tau_h3 = zeros(1,100000);
n_eq3 = zeros(1,100000);
m_eq3 = zeros(1,100000);
h_eq3 = zeros(1,100000);

tau_n4 = zeros(1,100000);
tau_m4 = zeros(1,100000);
tau_h4 = zeros(1,100000);
n_eq4 = zeros(1,100000);
m_eq4 = zeros(1,100000);
h_eq4 = zeros(1,100000);

V1(1) = V1_0; %initial given value
V2(1) = V2_0;
V3(1) = V3_0;
V4(1) = V4_0;

n1(1) = 0.3;
m1(1) = 0.05;
h1(1) = 0.6;

n2(1) = 0.3;
m2(1) = 0.05;
h2(1) = 0.6;

n3(1) = 0.3;
m3(1) = 0.05;
h3(1) = 0.6;

n4(1) = 0.3;
m4(1) = 0.05;
h4(1) = 0.6;

%% set alphot & betot & tau-im & einsofim due to given equations %%

alpha_n1(1) = (0.01*(V1(1) + 55))/(1 - exp(-0.1*(V1(1) + 55)));
alpha_m1(1) = (0.1*(V1(1) + 40))/(1-exp(-0.1*(V1(1) + 40)));
alpha_h1(1) = 0.07*exp(-0.05*(V1(1) + 65));
beta_n1(1) = 0.125*exp(-0.0125*(V1(1) + 65));
beta_m1(1) = 4*exp(-0.0556*(V1(1) + 65));
beta_h1(1) = 1/(1 + exp(-0.1*(V1(1) + 35)));

alpha_n2(1) = (0.01*(V2(1) + 55))/(1 - exp(-0.1*(V2(1) + 55)));
alpha_m2(1) = (0.1*(V2(1) + 40))/(1-exp(-0.1*(V2(1) + 40)));
alpha_h2(1) = 0.07*exp(-0.05*(V2(1) + 65));
beta_n2(1) = 0.125*exp(-0.0125*(V2(1) + 65));
beta_m2(1) = 4*exp(-0.0556*(V2(1) + 65));
beta_h2(1) = 1/(1 + exp(-0.1*(V2(1) + 35)));

alpha_n3(1) = (0.01*(V3(1) + 55))/(1 - exp(-0.1*(V3(1) + 55)));
alpha_m3(1) = (0.1*(V3(1) + 40))/(1-exp(-0.1*(V3(1) + 40)));
alpha_h3(1) = 0.07*exp(-0.05*(V3(1) + 65));
beta_n3(1) = 0.125*exp(-0.0125*(V3(1) + 65));
beta_m3(1) = 4*exp(-0.0556*(V3(1) + 65));
beta_h3(1) = 1/(1 + exp(-0.1*(V3(1) + 35)));

alpha_n4(1) = (0.01*(V4(1) + 55))/(1 - exp(-0.1*(V4(1) + 55)));
alpha_m4(1) = (0.1*(V4(1) + 40))/(1-exp(-0.1*(V4(1) + 40)));
alpha_h4(1) = 0.07*exp(-0.05*(V4(1) + 65));
beta_n4(1) = 0.125*exp(-0.0125*(V4(1) + 65));
beta_m4(1) = 4*exp(-0.0556*(V4(1) + 65));
beta_h4(1) = 1/(1 + exp(-0.1*(V4(1) + 35)));

tau_n1(1) = 1/(alpha_n1(1) + beta_n1(1));
tau_m1(1) = 1/(alpha_m1(1) + beta_m1(1));
tau_h1(1) = 1/(alpha_h1(1) + beta_h1(1));
n_eq1(1) = (alpha_n1(1))/(alpha_n1(1) + beta_n1(1));
m_eq1(1) = (alpha_m1(1))/(alpha_m1(1) + beta_m1(1));
h_eq1(1) = (alpha_h1(1))/(alpha_h1(1) + beta_h1(1));

tau_n2(1) = 1/(alpha_n2(1) + beta_n2(1));
tau_m2(1) = 1/(alpha_m2(1) + beta_m2(1));
tau_h2(1) = 1/(alpha_h2(1) + beta_h2(1));
n_eq2(1) = (alpha_n2(1))/(alpha_n2(1) + beta_n2(1));
m_eq2(1) = (alpha_m2(1))/(alpha_m2(1) + beta_m2(1));
h_eq2(1) = (alpha_h2(1))/(alpha_h2(1) + beta_h2(1));

tau_n3(1) = 1/(alpha_n3(1) + beta_n3(1));
tau_m3(1) = 1/(alpha_m3(1) + beta_m3(1));
tau_h3(1) = 1/(alpha_h3(1) + beta_h3(1));
n_eq3(1) = (alpha_n3(1))/(alpha_n3(1) + beta_n3(1));
m_eq3(1) = (alpha_m3(1))/(alpha_m3(1) + beta_m3(1));
h_eq3(1) = (alpha_h3(1))/(alpha_h3(1) + beta_h3(1));

tau_n4(1) = 1/(alpha_n4(1) + beta_n4(1));
tau_m4(1) = 1/(alpha_m4(1) + beta_m4(1));
tau_h4(1) = 1/(alpha_h4(1) + beta_h4(1));
n_eq4(1) = (alpha_n4(1))/(alpha_n4(1) + beta_n4(1));
m_eq4(1) = (alpha_m4(1))/(alpha_m4(1) + beta_m4(1));
h_eq4(1) = (alpha_h4(1))/(alpha_h4(1) + beta_h4(1));

t = linspace(0,1000,100000); %time vector

%% set current %%

Ie_1 = zeros(1,100000);
Ie_2 = zeros(1,100000);
Ie_3 = zeros(1,100000);
Ie_4 = zeros(1,100000);

%% sim 1 %%

%Ie_1(10000:11000) = 3.63; 

%% sim 2 %%

Ie_1(10000:11000) = 2.81; 
Ie_2(10000:11000) = 2.81; 

%% calculate molichut %%

g_31 = (radius3*(radius1^2))/((r_L*L*(L*(radius1^2))) + (L*(radius3^2)));
g_32 = (radius3*(radius2^2))/((r_L*L*(L*(radius2^2))) + (L*(radius3^2)));
g_43 = (radius4*(radius3^2))/((r_L*L*(L*(radius3^2))) + (L*(radius4^2)));
g_13 = (radius1*(radius3^2))/((r_L*L*(L*(radius3^2))) + (L*(radius1^2)));
g_23 = (radius2*(radius3^2))/((r_L*L*(L*(radius3^2))) + (L*(radius2^2)));
g_34 = (radius3*(radius4^2))/((r_L*L*(L*(radius4^2))) + (L*(radius3^2)));

%% set mekadmim %%

A1 = zeros(1,100000); %tamid
A2 = zeros(1,100000); %tamid
A31 = (1/C)*g_31*ones(1,100000);
A32 = (1/C)*g_32*ones(1,100000);
A4 = (1/C)*g_43*ones(1,100000);

B1 = zeros(1,100000);
B2 = zeros(1,100000);
B31 = zeros(1,100000);
B32 = zeros(1,100000);
B4 = zeros(1,100000);

B1(1) = (-1/C)*(Gna*(m1(1)^3)*h1(1) + Gk*(n1(1)^4) + Gl + 0 + g_13);
B2(1) = (-1/C)*(Gna*(m2(1)^3)*h2(1) + Gk*(n2(1)^4) + Gl + 0 + g_23);
B31(1) = (-1/C)*(Gna*(m3(1)^3)*h3(1) + Gk*(n3(1)^4) + Gl + g_31 + g_34);
B32(1) = (-1/C)*(Gna*(m3(1)^3)*h3(1) + Gk*(n3(1)^4) + Gl + g_32 + g_34);
B4(1) = (-1/C)*(Gna*(m4(1)^3)*h4(1) + Gk*(n4(1)^4) + Gl + g_43 + 0);

C1 = (1/C)*g_13*ones(1,100000);
C2 = (1/C)*g_23*ones(1,100000);
C3 = (1/C)*g_34*ones(1,100000);
C4 = zeros(1,100000); %tamid

D1 = zeros(1,100000);
D2 = zeros(1,100000);
D3 = zeros(1,100000);
D4 = zeros(1,100000);

D1(1) = (1/C)*(Gna*(m1(1)^3)*h1(1)*Ena + Gk*(n1(1)^4)*Ek + Gl*El + Ie_1(1)/(2*pi*radius1*L));
D2(1) = (1/C)*(Gna*(m2(1)^3)*h2(1)*Ena + Gk*(n2(1)^4)*Ek + Gl*El + Ie_2(1)/(2*pi*radius2*L));
D3(1) = (1/C)*(Gna*(m3(1)^3)*h3(1)*Ena + Gk*(n3(1)^4)*Ek + Gl*El + Ie_3(1)/(2*pi*radius3*L));
D4(1) = (1/C)*(Gna*(m4(1)^3)*h4(1)*Ena + Gk*(n4(1)^4)*Ek + Gl*El + Ie_4(1)/(2*pi*radius4*L));

%% set mishtanei ezer from the book

a1 = zeros(1,100000); %tamid
a2 = zeros(1,100000); %tamid
a31 = dt*A31;
a32 = dt*A32;
a4 = dt*A4; 

b1 = zeros(1,100000);
b2 = zeros(1,100000);
b31 = zeros(1,100000);
b32 = zeros(1,100000);
b4 = zeros(1,100000);
b3 = zeros(1,100000);

c1 = dt*C1;
c2 = dt*C2;
c3 = dt*C3;
c4 = zeros(1,100000); %tamid

d1 = zeros(1,100000);
d2 = zeros(1,100000);
d31 = zeros(1,100000);
d32 = zeros(1,100000);
d4 = zeros(1,100000);
d3 = zeros(1,100000);

%% numerical solution RK 2nd order %%
 
for i = 2:100000
    %% set RK 2nd order for m,h,n and update them then for membrane from HW1 %%
    %% cell 1
    k1_m1 = dt*((m_eq1(i-1) - m1(i-1))/tau_m1(i-1));
    k2_m1 = dt*((m_eq1(i-1) - (m1(i-1) + k1_m1))/tau_m1(i-1));
    m1(i) = m1(i-1) + 0.5*(k1_m1 + k2_m1); %RK2 update rule
    
    k1_n1 = dt*((n_eq1(i-1) - n1(i-1))/tau_n1(i-1));
    k2_n1 = dt*((n_eq1(i-1) - (n1(i-1) + k1_n1))/tau_n1(i-1));
    n1(i) = n1(i-1) + 0.5*(k1_n1 + k2_n1); %RK2 update rule
    
    k1_h1 = dt*((h_eq1(i-1) - h1(i-1))/tau_h1(i-1));
    k2_h1 = dt*((h_eq1(i-1) - (h1(i-1) + k1_h1))/tau_h1(i-1));
    h1(i) = h1(i-1) + 0.5*(k1_h1 + k2_h1); %RK2 update rule
    
    %% cell 2
    k1_m2 = dt*((m_eq2(i-1) - m2(i-1))/tau_m2(i-1));
    k2_m2 = dt*((m_eq2(i-1) - (m2(i-1) + k1_m2))/tau_m2(i-1));
    m2(i) = m2(i-1) + 0.5*(k1_m2 + k2_m2); %RK2 update rule
    
    k1_n2 = dt*((n_eq2(i-1) - n2(i-1))/tau_n2(i-1));
    k2_n2 = dt*((n_eq2(i-1) - (n2(i-1) + k1_n2))/tau_n2(i-1));
    n2(i) = n2(i-1) + 0.5*(k1_n2 + k2_n2); %RK2 update rule
    
    k1_h2 = dt*((h_eq2(i-1) - h2(i-1))/tau_h2(i-1));
    k2_h2 = dt*((h_eq2(i-1) - (h2(i-1) + k1_h2))/tau_h2(i-1));
    h2(i) = h2(i-1) + 0.5*(k1_h2 + k2_h2); %RK2 update rule
   
    %% cell 3
    k1_m3 = dt*((m_eq3(i-1) - m3(i-1))/tau_m3(i-1));
    k2_m3 = dt*((m_eq3(i-1) - (m3(i-1) + k1_m3))/tau_m3(i-1));
    m3(i) = m3(i-1) + 0.5*(k1_m3 + k2_m3); %RK2 update rule
    
    k1_n3 = dt*((n_eq3(i-1) - n3(i-1))/tau_n3(i-1));
    k2_n3 = dt*((n_eq3(i-1) - (n3(i-1) + k1_n3))/tau_n3(i-1));
    n3(i) = n3(i-1) + 0.5*(k1_n3 + k2_n3); %RK2 update rule
    
    k1_h3 = dt*((h_eq3(i-1) - h3(i-1))/tau_h3(i-1));
    k2_h3 = dt*((h_eq3(i-1) - (h3(i-1) + k1_h3))/tau_h3(i-1));
    h3(i) = h3(i-1) + 0.5*(k1_h3 + k2_h3); %RK2 update rule
    
    %% cell 4
    k1_m4 = dt*((m_eq4(i-1) - m4(i-1))/tau_m4(i-1));
    k2_m4 = dt*((m_eq4(i-1) - (m4(i-1) + k1_m4))/tau_m4(i-1));
    m4(i) = m4(i-1) + 0.5*(k1_m4 + k2_m4); %RK2 update rule
    
    k1_n4 = dt*((n_eq4(i-1) - n4(i-1))/tau_n4(i-1));
    k2_n4 = dt*((n_eq4(i-1) - (n4(i-1) + k1_n4))/tau_n4(i-1));
    n4(i) = n4(i-1) + 0.5*(k1_n4 + k2_n4); %RK2 update rule
    
    k1_h4 = dt*((h_eq4(i-1) - h4(i-1))/tau_h4(i-1));
    k2_h4 = dt*((h_eq4(i-1) - (h4(i-1) + k1_h4))/tau_h4(i-1));
    h4(i) = h4(i-1) + 0.5*(k1_h4 + k2_h4); %RK2 update rule
    
    %% update current mekadmim 
    
    B1(i) = (-1/C)*(Gna*(m1(i)^3)*h1(i) + Gk*(n1(i)^4) + Gl + 0 + g_13);
    B2(i) = (-1/C)*(Gna*(m2(i)^3)*h2(i) + Gk*(n2(i)^4) + Gl + 0 + g_23);
    B31(i) = (-1/C)*(Gna*(m3(i)^3)*h3(i) + Gk*(n3(i)^4) + Gl + g_31 + g_34);
    B32(i) = (-1/C)*(Gna*(m3(i)^3)*h3(i) + Gk*(n3(i)^4) + Gl + g_32 + g_34);
    B4(i) = (-1/C)*(Gna*(m4(i)^3)*h4(i) + Gk*(n4(i)^4) + Gl + g_43 + 0);
    
    D1(i) = (1/C)*(Gna*(m1(i)^3)*h1(i)*Ena + Gk*(n1(i)^4)*Ek + Gl*El + Ie_1(i)/(2*pi*radius1*L));
    D2(i) = (1/C)*(Gna*(m2(i)^3)*h2(i)*Ena + Gk*(n2(i)^4)*Ek + Gl*El + Ie_2(i)/(2*pi*radius2*L));
    D3(i) = (1/C)*(Gna*(m3(i)^3)*h3(i)*Ena + Gk*(n3(i)^4)*Ek + Gl*El + Ie_3(i)/(2*pi*radius3*L));
    D4(i) = (1/C)*(Gna*(m4(i)^3)*h4(i)*Ena + Gk*(n4(i)^4)*Ek + Gl*El + Ie_4(i)/(2*pi*radius4*L));
    
    b1(i) = dt*B1(i);
    b2(i) = dt*B2(i);
    b31(i) = dt*B31(i);
    b32(i) = dt*B32(i);
    b4(i) = dt*B4(i);
    b3(i) = 0.5*(b31(i) + b32(i)); %normelize their affection
    
    d1(i) = (D1(i) + 0 + B1(i)*V1(i-1) + C1(i)*V3(i-1))*dt;
    d2(i) = (D2(i) + 0 + B2(i)*V2(i-1) + C2(i)*V3(i-1))*dt;
    d31(i) = (D3(i) + A31(i)*V1(i-1) + B31(i)*V3(i-1) + C3(i)*V4(i-1))*dt;
    d32(i) = (D3(i) + A32(i)*V2(i-1) + B32(i)*V3(i-1) + C3(i)*V4(i-1))*dt;
    d4(i) = (D4(i) + A4(i)*V3(i-1) + B4(i)*V4(i-1) + 0)*dt;
    d3(i) = 0.5*(d31(i) + d32(i)); %normelize their affection
    
    %% update deltot voltage according to Theoritical neuroscience book p 225-227
    
    delta_V41 = (d4(i) + (a4(i)*(d31(i)+(a31(i)*d31(i))/(1-b1(i))))/...
        (1 - (b31(i) + (a31(i)*c1(i))/(1 - b1(i)))))/...
        (1 - (b4(i) + (a4(i)*c3(i))/...
        (1 - (b31(i) + (a31(i)*c1(i))/(1 - b1(i))))));
    delta_V42 = (d4(i) + (a4(i)*(d32(i)+(a32(i)*d32(i))/(1-b2(i))))/...
        (1 - (b32(i) + (a32(i)*c2(i))/(1 - b2(i)))))/...
        (1 - (b4(i) + (a4(i)*c3(i))/...
        (1 - (b32(i) + (a32(i)*c2(i))/(1 - b2(i))))));
    delta_V4 = 0.5*(delta_V41 + delta_V42);
    
    delta_V32 = (c3(i)*delta_V4 + d32(i) + (a32(i)*(d2(i)))/(1 - b2(i)))/...
        (1 - (b32(i) + (a32(i)*c2(i))/(1 - b2(i))));
    delta_V31 = (c3(i)*delta_V4 + d31(i) + (a31(i)*(d1(i)))/(1 - b1(i)))/...
        (1 - (b31(i) + (a31(i)*c1(i))/(1 - b1(i))));
        
    delta_V2 = (c2(i)*delta_V32 + d2(i))/(1 - b2(i));
    delta_V1 = (c1(i)*delta_V31 + d1(i))/(1 - b1(i));
   
    %% update voltage according to Theoritical neuroscience book p 225-227
    
    V1(i) = V1(i-1) + delta_V1;
    V2(i) = V2(i-1) + delta_V2;
    V3(i) = V3(i-1) + 0.5*(delta_V31 + delta_V32);
    V4(i) = V4(i-1) + delta_V4;
    
    %%% update next alphot & betot & tau-im & einsofim %%%
    alpha_n1(i) = (0.01*(V1(i) + 55))/(1 - exp(-0.1*(V1(i) + 55)));
    alpha_m1(i) = (0.1*(V1(i) + 40))/(1-exp(-0.1*(V1(i) + 40)));
    alpha_h1(i) = 0.07*exp(-0.05*(V1(i) + 65));
    beta_n1(i) = 0.125*exp(-0.0125*(V1(i) + 65));
    beta_m1(i) = 4*exp(-0.0556*(V1(i) + 65));
    beta_h1(i) = 1/(1 + exp(-0.1*(V1(i) + 35)));
    
    alpha_n2(i) = (0.01*(V2(i) + 55))/(1 - exp(-0.1*(V2(i) + 55)));
    alpha_m2(i) = (0.1*(V2(i) + 40))/(1-exp(-0.1*(V2(i) + 40)));
    alpha_h2(i) = 0.07*exp(-0.05*(V2(i) + 65));
    beta_n2(i) = 0.125*exp(-0.0125*(V2(i) + 65));
    beta_m2(i) = 4*exp(-0.0556*(V2(i) + 65));
    beta_h2(i) = 1/(1 + exp(-0.1*(V2(i) + 35)));
    
    alpha_n3(i) = (0.01*(V3(i) + 55))/(1 - exp(-0.1*(V3(i) + 55)));
    alpha_m3(i) = (0.1*(V3(i) + 40))/(1-exp(-0.1*(V3(i) + 40)));
    alpha_h3(i) = 0.07*exp(-0.05*(V3(i) + 65));
    beta_n3(i) = 0.125*exp(-0.0125*(V3(i) + 65));
    beta_m3(i) = 4*exp(-0.0556*(V3(i) + 65));
    beta_h3(i) = 1/(1 + exp(-0.1*(V3(i) + 35)));
    
    alpha_n4(i) = (0.01*(V4(i) + 55))/(1 - exp(-0.1*(V4(i) + 55)));
    alpha_m4(i) = (0.1*(V4(i) + 40))/(1-exp(-0.1*(V4(i) + 40)));
    alpha_h4(i) = 0.07*exp(-0.05*(V4(i) + 65));
    beta_n4(i) = 0.125*exp(-0.0125*(V4(i) + 65));
    beta_m4(i) = 4*exp(-0.0556*(V4(i) + 65));
    beta_h4(i) = 1/(1 + exp(-0.1*(V4(i) + 35)));
    
    tau_n1(i) = 1/(alpha_n1(i) + beta_n1(i));
    tau_m1(i) = 1/(alpha_m1(i) + beta_m1(i));
    tau_h1(i) = 1/(alpha_h1(i) + beta_h1(i));
    n_eq1(i) = (alpha_n1(i))/(alpha_n1(i) + beta_n1(i));
    m_eq1(i) = (alpha_m1(i))/(alpha_m1(i) + beta_m1(i));
    h_eq1(i) = (alpha_h1(i))/(alpha_h1(i) + beta_h1(i));
    
    tau_n2(i) = 1/(alpha_n2(i) + beta_n2(i));
    tau_m2(i) = 1/(alpha_m2(i) + beta_m2(i));
    tau_h2(i) = 1/(alpha_h2(i) + beta_h2(i));
    n_eq2(i) = (alpha_n2(i))/(alpha_n2(i) + beta_n2(i));
    m_eq2(i) = (alpha_m2(i))/(alpha_m2(i) + beta_m2(i));
    h_eq2(i) = (alpha_h2(i))/(alpha_h2(i) + beta_h2(i));
    
    tau_n3(i) = 1/(alpha_n3(i) + beta_n3(i));
    tau_m3(i) = 1/(alpha_m3(i) + beta_m3(i));
    tau_h3(i) = 1/(alpha_h3(i) + beta_h3(i));
    n_eq3(i) = (alpha_n3(i))/(alpha_n3(i) + beta_n3(i));
    m_eq3(i) = (alpha_m3(i))/(alpha_m3(i) + beta_m3(i));
    h_eq3(i) = (alpha_h3(i))/(alpha_h3(i) + beta_h3(i));
    
    tau_n4(i) = 1/(alpha_n4(i) + beta_n4(i));
    tau_m4(i) = 1/(alpha_m4(i) + beta_m4(i));
    tau_h4(i) = 1/(alpha_h4(i) + beta_h4(i));
    n_eq4(i) = (alpha_n4(i))/(alpha_n4(i) + beta_n4(i));
    m_eq4(i) = (alpha_m4(i))/(alpha_m4(i) + beta_m4(i));
    h_eq4(i) = (alpha_h4(i))/(alpha_h4(i) + beta_h4(i));
    
end  
    
figure;
subplot(2,2,1); plot(t, V1); title('V1 vs time');
xlabel('t [msec]'); ylabel('V1 [mV]'); 
subplot(2,2,2); plot(t, V2); title('V2 vs time');
xlabel('t [msec]'); ylabel('V2 [mV]');  
subplot(2,2,3); plot(t, V3); title('V3 vs time');
xlabel('t [msec]'); ylabel('V3 [mV]'); 
subplot(2,2,4); plot(t, V4); title('V4 vs time');
xlabel('t [msec]'); ylabel('V4 [mV]'); 
suptitle(strcat('I1 =', num2str(Ie_1(10500)),' [nA]', '   I2 =', num2str(Ie_2(10500)),' [nA]'));