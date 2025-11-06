function [Y,sim_Ftotal_cycle,integfail_flag,cycle_fail_index,complex_flag]=cycle_1_pi_H(params,iemg_profile,Pi_set,H_set,met_dyn_set)
    % Experimental data
    data_resting=readtable("../raw_data/initial_state.xlsx"); % resting levels of state variable
    % Set metabolite concentrations, 
    MgATP = 8.2; % Assumption made in Broxtermann et al 2017  
    MgADP = data_resting{1,2}*10^-3; 
    %Pi = data_resting{4,2}; % time 0 value reported in Broxtermann et al 2017
    Pi = Pi_set;
    %Pi = 0; % Experimentally estimated resting levels by Umass team
    Pcr = data_resting{2,2};% Experimentally estimated resting levels by Umass team
    %Pcr = 40; % Experimentally estimated resting levels by Umass team
    SL0 = 3.23;%2.2; % um [114 mm = 1.3758 um; 116 mm = 1.5869; 118mm = 1.8114; 120mm = 2.0403; 128 mm = 2.8346; 130mm = 2.9728; 132mm= 3.0980]
    H = H_set;
    %H = data_resting{3,2}; % Experimentally estimated resting levels by Umass team
    init = [zeros(1,9),SL0, Pi,MgADP, Pcr,H,MgATP]; % Initial conditions for the model
    cycle_time=3;% 3s contraction 2s relaxation in Broxtermann et al., 2017; 
    tspan = 0:0.1:cycle_time;
    n=length(tspan);
    integfail_flag=0;
    cycle_fail_index=0;
    complex_flag = 0;
    SL_set=3.23;
    iemg = iemg_profile;
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
    [T, Y] = ode15s(@Model_XB_human_QC_metdyn_set,tspan,init,options,SL_set,params,iemg,Pcr,met_dyn_set);
    sim_Ftotal_cycle = zeros(n,1);
    for j=1:n
        [~, sim_Ftotal_cycle(j),~,~,~,~,~,~,~,~,~,~,~,~,~] = Model_XB_human_QC(T(j),Y(j,:),SL_set,params,iemg,Pcr);
    end
    k=length(T);
    if k<n
        integfail_flag=1; 
    end
    if ~isreal(Y)
        complex_flag =1;
    end 
end