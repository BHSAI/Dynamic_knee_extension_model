function [ dYdT, Ftotal, F_active,F_passive,F_active_s,sov_thick, B_process, C_process, dCK, dK3, dCK_r, dPi_cons, dH_cons, dGly, dAdk] = Model_XB_human_QC(t,y,SLset,par,Pu0,Pcr0)
    
    alpha1 = 1*par(8); % 1/um
    alpha2 = 1*par(9); % 1/um
    alpha3 = 1*par(10);% 1/um
    s3 = 1*par(11);      % um
    %K_Pi = 3.0*par(12);    % mM Pi dissociation constant
    %K_Pi = 3.0*7.07020430862858;    % mM Pi dissociation constant estimated by fitting Broxterman et al 2017
    K_Pi = par(19);
    %K_T = 1*par(15);     % mM MgATP dissociation constant
    K_T = par(20);     % mM MgATP dissociation constant estimated by fitting Broxterman et al 2017 
    %K_D = 0.4*0.194; % MgADP dissociation constant from Yamashita etal (Circ Res. 1994; 74:1027-33).
    % K_D = 1*0.07; % mM
    %K_D = par(20); % mM
    K_D = par(21); % mM MgADP dissociation constant estimated by fitting Broxterman et al 2017
    %eta = par(24);
    kf1 = par(14); % rate of PCr consumption
    kf2 = par(15); % rate of Pi consumption :without h1 binding polynomial, 0.004
    %Pcr0 = 6; 
    %kad_d = 0.12; % adp consumption
    %kad_d=par(21);
    %Kh1 = 3e4*10^-7; % mM without h1 binding polynomial, 6e4*10^-7
    %Kh1 = par(18)*10^-7;
    Kh1 = par(22); %mM proton dissociation constant estimated by fitting Broxterman et al 2017  
    %kf3 = 7.2e1; % rate of H+ consumption 1/sec
    kf_1 = par(16);
    k_gly = par(17);
    k_adk = par(18);
    dSL_set = par(23);
    % H = 10^-pH;
    %% State Variables
    P1o = y(1);
    P1i = y(2);
    P1w = y(3);
    P2o = y(4);
    P2i = y(5);
    P2w = y(6);
    P3o = y(7);
    P3i = y(8);
    P3w = y(9);
    SL = y(10);
    Pi = y(11);
    MgADP = y(12);
    Pcr = y(13);
    H = y(14); % 10^-pH
    MgATP=y(15);
    Pu = Pu0 - P1o - P2o - P3o;
    % K_D = K_D*(1/(1+ H/Kh1)); %% Simple modification of ADP dissociation
    % based on protons
    g1 = (MgADP/K_D)/(1 + MgADP/K_D + MgATP/K_T); 
    %g1 = (MgADP/K_D)/(1 + MgADP/K_D + MgATP/K_T + (MgADP/K_D)*H/Kh1); 
    g2 = (MgATP/K_T)/(1 + MgATP/K_T + MgADP/K_D ); 
    %g2 = (MgATP/K_T)/(1 + MgATP/K_T + MgADP/K_D + (MgADP/K_D)*H/Kh1); 
    f1 = (Pi/K_Pi)/(1 + Pi/K_Pi); f2 = 1/(1 + Pi/K_Pi); 
    %f1 = (Pi/K_Pi)/(1 + Pi/K_Pi + H/Kh1); f2 = 1/(1 + Pi/K_Pi + H/Kh1); 
    % f1 = (Pi/K_Pi)/(1 + Pi/K_Pi + (Pi/K_Pi)*(H/Kh1)); f2 = 1/(1 + Pi/K_Pi + (Pi/K_Pi)*(H/Kh1)); 
    %h1 = 1/(1 + H/Kh1);    
    %h1=1;
    h1=(H/Kh1)/(1 + H/Kh1); h2 = 1/(1 + H/Kh1);

    % rate parameter                          % Units

    kf = par(1);      % 1/sec
    kb = par(2)*f1;   % 1/sec
    k1 = par(3)*f2;   % 1/sec
    %k_1 =rho_myosin*par(4)*Q10s(1)^((TmpC-17)/10);     % 1/sec
    k_1 =par(4)*h1;     % 1/sec
    %k2 = rho_myosin*par(5)*1*Q10s(2)^((TmpC-17)/10);    % 1/sec
    k2 = par(5)*h2;    % 1/sec
    k_2 = par(6)*1*g1;% 1/sec
    k3 = par(7)*g2;   % 1/sec
    % kpe1 = par(13)*Q10s(3)^((TmpC-17)/10);
    %eta =1e-3*par(14)*Q10s(3)^((TmpC-17)/10);   % N.sec/(mm^2.um)
    kstiff1 = 1*par(12); % mN/(mm^2.um)
    % kpe2 = par(16)*Q10s(5)^((TmpC-17)/10);  
    kstiff2 = 1*par(13); % mN/(mm^2.um)
    
    %% Sarcomere geometry (um)
    SL_max = 4.2; %2.4
    SL_min = 1.4;
    SL_rest = 1.9;  % (um)
    L_thick = 1.65;
    L_hbare = 0.12;
    % L_thin = 1.2;
    L_thin = 1.27;
    %% Thin filament activation rates
    % Sarcomere geometry
    sovr_ze = min(L_thick*0.50, SL*0.38);
    % sovr_cle = max(SL*0.5 - (SL-L_thin),L_hbare*0.5);
    sovr_cle = max(abs(SL*0.50 - (SL-L_thin)),L_hbare*0.4);
    L_sovr = sovr_ze - sovr_cle; % Length of single overlap region
        
    % Overlap fraction for thick filament
    sov_thick = L_sovr*2/(L_thick - L_hbare);
    % Overlap fraction for thin filament
    % sov_thin = L_sovr/L_thin;
    % figure; plot(SL,sov_thick,'r')
    %% Stretch-sensitive rates
    f_alpha1o = (P1o - alpha1*P1i + 0.5*(alpha1*alpha1)*P1w);
    f_alpha1i = (P1i - alpha1*P1w);
    
    alpha0 = 1*alpha1;
    f_alpha0o = (P2o + alpha0*P2i + 0.5*alpha0*alpha0*P2w);
    f_alpha0i = (P2i + alpha0*P2w);
    
    f_alpha2o = (P2o - alpha2*P2i + 0.5*(alpha2*alpha2)*P2w);
    f_alpha2i = (P2i - alpha2*P2w);
    
    alpha2b = 0; 
    f_alphao = (P3o + alpha2b*P3i + 0.5*(alpha2b*alpha2b)*P3w);
    f_alphai = (P3i + alpha2b*P3w);
    
    f_alpha3o = (P3o + alpha3*(s3*s3*P3o + 2*s3*P3i + P3w));
    f_alpha3i = (P3i + alpha3*(s3*s3*P3i + 2*s3*P3w));
    
    %% Compute Active & Passive Force
    % Active Force
    dr = 0.01; % Power-stroke Size; Units: um
    B_process = 0.6*kstiff2*dr*P3o;   % Force due to XB cycling (without h1 it is 0.56
    % B_process = 0.46*kstiff2*dr*P3o;   % Use it for simultion K_D; Force due to XB cycling (without h1 it is 0.56
    C_process = kstiff1*(P2i+P3i);% Force due to stretching of XBs
    F_active_s = sov_thick*(B_process + C_process); % mN/mm^2
    % convert to units of mN/mm^2 using the cross sectional area (CSA)
    % rabbit GM csa 8.63 cm^2
    % CSA = 10*8.63e-2; % m^2
    % physiological cross sectional area
    % Soleus PCSA = 98 (sd=9) cm^2 (Sopher et al JSA 2017)(% ranges from
    % 80 to 230 cm^2 from different sources)
    % CSA = 118*1e2; % mm^2 soleus
    CSA = 100*1e2; % mm^2 Quadriceps
    F_active = 1e-3*F_active_s*CSA; % N
    
    % Non-linear Passive force; Adopted from Rice etal (Biophys J. 2008 Sep;95(5):2368-90)
    [F_preload,~] = passiveForces_rabbit_QC(SLset,SLset);
    [F_passive,FSEE] = passiveForces_rabbit_QC(SL,SLset);
    
    f_myofibril = 0.60; % Percent Myofibril in Skeletal Muscle from Haun etal [PMID: 30930796]
    Ftotal = f_myofibril*(F_active +  F_passive);%
    
    %% XB ODEs
    % dSL = ((intf/eta) - dfxb/kpe1)*heav(SL-SL_min)*heav(SL_max-SL)/den;
    %eta_n = eta*CSA; % N*s/um; 
    %intf = (-Ftotal + F_preload + FSEE); 
    %dSL = (intf/eta_n)*heav(SL-SL_min)*heav(SL_max-SL); 
    dSL = dSL_set; 
    dP1o = kf*Pu   - kb*P1o - k1*f_alpha1o + k_1*f_alpha0o;
    dP1i = 1*dSL*P1o - kb*P1i - k1*f_alpha1i + k_1*f_alpha0i;
    dP1w = 2*dSL*P1i - kb*P1w - k1*P1w + k_1*P2w;
    
    dP2o =         - k_1*f_alpha0o - k2*f_alpha2o + k_2*f_alphao + k1*f_alpha1o;
    dP2i = 1*dSL*P2o - k_1*f_alpha0i - k2*f_alpha2i + k_2*f_alphai + k1*f_alpha1i;
    dP2w = 2*dSL*P2i - k_1*P2w       - k2*P2w + k_2*P3w + k1*P1w;
    
    dP3o =         + k2*f_alpha2o - k_2*f_alphao - k3*f_alpha3o;
    dP3i = 1*dSL*P3o + k2*f_alpha2i - k_2*f_alphai - k3*f_alpha3i;
    dP3w = 2*dSL*P3i + k2*P2w       - k_2*P3w - k3*P3w;
    
    % metabolic network
    dPi =  650*1e-3*k3*f_alpha3o - 1.0*kf2*Pi - k_gly*MgADP*Pi;
    dPi_cons=1.0*kf2*Pi;
    gamma = 0.6; %stoichimetry of proton release from ATP hydrolysis taken from Barclay 2017
    Ka=1000*10^(-6.75);beta = (2.3*Ka*4.7*H)/((Ka+H)^2); delta_pH=7-(-log10((10^-3)*H));
    dH = gamma*650*1e-3*k3*f_alpha3o -kf1*Pcr*MgADP- beta*delta_pH + k_gly*MgADP*Pi + kf_1*(Pcr0-Pcr)*MgATP;
    dH_cons= beta*delta_pH;
    dMgADP = 650*1e-3*k3*f_alpha3o -kf1*Pcr*MgADP -k_gly*MgADP*Pi + kf_1*(Pcr0-Pcr)*MgATP - k_adk*MgADP*MgADP;
    dCK=kf1*Pcr*MgADP;
    dK3=650*1e-3*k3*f_alpha3o;
    dGly=k_gly*MgADP*Pi; %glycolysis
    dAdk=k_adk*MgADP*MgADP; %adenylate kinase
    dPcr = -kf1*Pcr*MgADP + kf_1*(Pcr0-Pcr)*MgATP;
    dCK_r = kf_1*(Pcr0-Pcr)*MgATP;
    dMgATP=-dMgADP;
    dYdT = [dP1o; dP1i; dP1w; dP2o; dP2i; dP2w; dP3o; dP3i; dP3w; dSL; dPi; dMgADP; dPcr; dH; dMgATP];
end 