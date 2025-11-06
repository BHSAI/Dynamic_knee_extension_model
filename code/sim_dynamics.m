function [error,S,sim_Ftotal_cycles,R,integfail_flag,cycle_fail_index,complex_flag]=sim_dynamics(params)
    % fit iemg to hill equation and get the hill parameters
    iemg_table=readtable("../raw_data/iEMG_data_burnley.xlsx");
    iemg_time = iemg_table{:,1}; iemg_val = iemg_table{:,2};
    [p_iemg,~,~]=poly2(iemg_time,iemg_val);
    % Experimental data
    data_Pcr  = readtable('../raw_data/PCr_data_broxterman.csv'); % Pcr
    data_Pi  = readtable('../raw_data/Pi_data_broxterman.csv'); % Pi 
    %data_ADP  = readtable('data_exp/ADP_data_broxterman.xlsx'); % ADP 
    data_proton  = readtable('../raw_data/proton_data_broxterman.xlsx'); % proton
    data_force = readtable('../raw_data/force_data_broxterman.csv'); % force
    data_resting=readtable("../raw_data/initial_state.xlsx"); % resting levels of state variables
    cycle_index_exp=1:60; % Broxtermann et al 2017 conducted 60 maximum voluntary contraction 
    cycle_start_time=1:5:300;
    p1 = data_Pcr{2:end,2};
    p2 = data_Pi{2:end,2};
    %p3 = 10^-3*data_ADP{2:end,2};
    p3 = data_proton{2:end,2};
    p4 = data_force{1:end,2};
    % Set metabolite concentrations, 
    MgATP = 8.2; % Assumption made in Broxtermann et al 2017  
    MgADP = data_resting{1,2}*10^-3; 
    Pi = data_resting{4,2}; % time 0 value reported in Broxtermann et al 2017 
    %Pi = 0; % Experimentally estimated resting levels by Umass team
    Pcr = data_resting{2,2};% Experimentally estimated resting levels by Umass team
    %Pcr = 40; % Experimentally estimated resting levels by Umass team
    SL0 = 3.23;%2.2; % um [114 mm = 1.3758 um; 116 mm = 1.5869; 118mm = 1.8114; 120mm = 2.0403; 128 mm = 2.8346; 130mm = 2.9728; 132mm= 3.0980]
    H = data_resting{3,2}; % Experimentally estimated resting levels by Umass team
    init = [zeros(1,9),SL0, Pi,MgADP, Pcr,H,MgATP]; % Initial conditions for the model
    cycles=1:1:max(cycle_index_exp);
    cycle_time=3;% 3s contraction 2s relaxation in Broxtermann et al., 2017; 
    tspan = 0:0.1:cycle_time;
    n=length(tspan);
    m=length(cycles);
    pi_p=zeros(m,1);
    Pu_init = zeros(m,3);
    Pu_fin = zeros(m,3);
    ADP_p=zeros(m,1);
    Pcr_p=zeros(m,1);
    H_p=zeros(m,1);
    ATP_p=zeros(m,1);
    sim_Ftotal_1s=zeros(m,1);
    sim_Ftotal_2s=zeros(m,1);
    sim_Ftotal_3s=zeros(m,1);
    sim_Ftotal = zeros(m,1);
    dCK_final = zeros(m,1);
    dK3_final = zeros(m,1);
    dCK_r_final = zeros(m,1);
    dPi_cons_final = zeros(m,1);
    dH_cons_final = zeros(m,1);
    dGly_cons_final = zeros(m,1);
    dAdk_cons_final = zeros(m,1); 
    sim_Ftotal_cycles=zeros(n,m);
    dCK_cycle=zeros(n,m);
    dK3_cycle=zeros(n,m);
    dCK_r_cycle=zeros(n,m);
    dPi_cons_cycle=zeros(n,m); 
    dH_cons_cycle=zeros(n,m); 
    dGly_cons_cycle = zeros(n,m); 
    dAdk_cons_cycle = zeros(n,m); 
    integfail_flag=0;
    cycle_fail_index=0;
    complex_flag = 0;
    for i=1:m
        SL_set=3.23;
        iemg=((p_iemg(1)*(cycle_start_time(i)^2))+(p_iemg(2)*(cycle_start_time(i)^1))+p_iemg(3))/100;
        options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
        [T, Y] = ode15s(@Model_XB_human_QC,tspan,init,options,SL_set,params,iemg,Pcr);
        k=length(T);
        if k<n
            integfail_flag=1; 
        end
        if integfail_flag==1
            cycle_fail_index=i;
        break   
        end
        init(11)=Y(n,11);%Pi
        init(12)=Y(n,12);%ADP
        init(13)=Y(n,13);%Pcr
        init(14)=Y(n,14);%H
        init(15)=Y(n,15);%ATP
        %Pu_init(i) =iemg-Y(1,1)+Y(1,2)+Y(1,6);
        %Pu_fin(i) =iemg-Y(n,1)+Y(n,2)+Y(n,6);
        Pu_init(i,:) =[Y(1,1) Y(1,4) Y(1,7)];
        Pu_fin(i,:) =[Y(n,1) Y(n,4) Y(n,7)];
        pi_p(i)=Y(n,11);
        ADP_p(i)=Y(n,12);
        Pcr_p(i)=Y(n,13);
        H_p(i)=Y(n,14);
        ATP_p(i)=Y(n,15);
        for j=1:n
            [~, sim_Ftotal_cycles(j,i),~,~,~,~,~,~,dCK_cycle(j,i),dK3_cycle(j,i),dCK_r_cycle(j,i),dPi_cons_cycle(j,i),dH_cons_cycle(j,i),dGly_cons_cycle(j,i),dAdk_cons_cycle(j,i)] = Model_XB_human_QC(T(j),Y(j,:),SL_set,params,iemg,Pcr);
        end
        sim_Ftotal(i) = max(sim_Ftotal_cycles(:,i));
        dCK_final(i) = dCK_cycle(n,i);
        dK3_final(i) = dK3_cycle(n,i);
        dCK_r_final(i) = dCK_r_cycle(n,i);
        dPi_cons_final(i) = dPi_cons_cycle(n,i);
        dH_cons_final(i) = dH_cons_cycle(n,1);
        dGly_cons_final(i) = dGly_cons_cycle(n,i);
        dAdk_cons_final(i) = dAdk_cons_cycle(n,i);
        [~, sim_Ftotal_3s(i),~,~,~,~,~,~,~,~,~,~,~,~,~] = Model_XB_human_QC(T(n),Y(n,:),SL_set,params,iemg,Pcr);
        [~, sim_Ftotal_2s(i),~,~,~,~,~,~,~,~,~,~,~,~,~] = Model_XB_human_QC(T(21),Y(21,:),SL_set,params,iemg,Pcr);
        [~, sim_Ftotal_1s(i),~,~,~,~,~,~,~,~,~,~,~,~,~] = Model_XB_human_QC(T(11),Y(11,:),SL_set,params,iemg,Pcr);
    end
    %S=[pi_p ADP_p Pcr_p H_p ATP_p sim_Ftotal];
    sim_Ftotal_2 = mean(horzcat(sim_Ftotal_1s,sim_Ftotal_2s,sim_Ftotal_3s),2);
    S=[pi_p ADP_p Pcr_p H_p ATP_p sim_Ftotal_2];
    R=[dCK_final dK3_final dCK_r_final dPi_cons_final dH_cons_final dGly_cons_final dAdk_cons_final];
    if ~isreal(S)
        complex_flag =1;
    end
    err1 = sum(((Pcr_p-p1)/max(p1)).^2);
    err2 = sum(((pi_p-p2)/max(p2)).^2);
    %err3 = sum(((ADP_p-p3)/max(p3)).^2);
    err3 = sum(((H_p-p3)/max(p3)).^2);
    err4 = sum(((sim_Ftotal_3s-p4)/max(p4)).^2);
    err5 = sum(((sim_Ftotal_2s-p4)/max(p4)).^2);
    err6 = sum(((sim_Ftotal_1s-p4)/max(p4)).^2);
    error = err1+err2+err3+err4+err5+err6;  
end