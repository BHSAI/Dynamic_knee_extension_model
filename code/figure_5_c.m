% read the param
param_table_dke=readtable('params/params.xlsx');
params= param_table_dke.Values;
% fit iemg to hill equation and get the hill parameters
iemg_table=readtable("../raw_data/iEMG_data_burnley.xlsx");
iemg_time = iemg_table{:,1}; iemg_val = iemg_table{:,2};
[p_iemg,~,~]=poly2(iemg_time,iemg_val);
% Experimental data
data_resting=readtable("../raw_data/initial_state.xlsx"); % resting levels of state variables
cycle_index_exp=1:60; % Broxtermann et al 2017 conducted 60 maximum voluntary contraction 
cycle_start_time=1:5:300;
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
n = length(tspan);
m = length(cycles);
pi_p = zeros(m,1);
ADP_p = zeros(m,1);
Pcr_p = zeros(m,1);
H_p = zeros(m,1);
ATP_p = zeros(m,1);
P1o_p = zeros(m,1);
P2o_p = zeros(m,1);
P2i_p = zeros(m,1);
P3o_p = zeros(m,1);
P3i_p = zeros(m,1);
sim_Ftotal_1s = zeros(m,1);
sim_Ftotal_2s = zeros(m,1);
sim_Ftotal_3s = zeros(m,1);
sim_F_passive_1s = zeros(m,1);
sim_F_passive_2s = zeros(m,1);
sim_F_passive_3s = zeros(m,1);
sim_sov_thick_1s = zeros(m,1);
sim_sov_thick_2s = zeros(m,1);
sim_sov_thick_3s = zeros(m,1);
sim_B_process_1s = zeros(m,1);
sim_B_process_2s = zeros(m,1);
sim_B_process_3s = zeros(m,1);
sim_C_process_1s = zeros(m,1);
sim_C_process_2s = zeros(m,1);
sim_C_process_3s = zeros(m,1);
%sim_Ftotal = zeros(m,1);
sim_Ftotal_cycles=zeros(n,m);
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
    P1o_p(i) = Y(n,1);
    P2o_p(i) = Y(n,4);
    P2i_p(i) = Y(n,5);
    P3o_p(i) = Y(n,7);
    P3i_p(i) = Y(n,8);
    pi_p(i)=Y(n,11);
    ADP_p(i)=Y(n,12);
    Pcr_p(i)=Y(n,13);
    H_p(i)=Y(n,14);
    ATP_p(i)=Y(n,15);
    for j=1:n
        [~, sim_Ftotal_cycles(j,i),~,~,~,~,~,~,~,~,~,~,~,~,~] = Model_XB_human_QC(T(j),Y(j,:),SL_set,params,iemg,Pcr);
    end
    %sim_Ftotal(i) = max(sim_Ftotal_cycles(:,i));
    [~, sim_Ftotal_3s(i),~,sim_F_passive_3s(i),~,sim_sov_thick_3s(i),sim_B_process_3s(i),sim_C_process_3s(i),~,~,~,~,~,~,~] = Model_XB_human_QC(T(n),Y(n,:),SL_set,params,iemg,Pcr);
    [~, sim_Ftotal_2s(i),~,sim_F_passive_2s(i),~,sim_sov_thick_2s(i),sim_B_process_2s(i),sim_C_process_2s(i),~,~,~,~,~,~,~] = Model_XB_human_QC(T(21),Y(21,:),SL_set,params,iemg,Pcr);
    [~, sim_Ftotal_1s(i),~,sim_F_passive_1s(i),~,sim_sov_thick_1s(i),sim_B_process_1s(i),sim_C_process_1s(i),~,~,~,~,~,~,~] = Model_XB_human_QC(T(11),Y(11,:),SL_set,params,iemg,Pcr);
end
%S=[pi_p ADP_p Pcr_p H_p ATP_p sim_Ftotal];
sim_Ftotal_2 = mean(horzcat(sim_Ftotal_1s,sim_Ftotal_2s,sim_Ftotal_3s),2);
sim_F_passive = mean(horzcat(sim_F_passive_1s,sim_F_passive_2s,sim_F_passive_3s),2);
sim_sov_thick = mean(horzcat(sim_sov_thick_1s,sim_sov_thick_2s,sim_sov_thick_3s),2);
sim_B_process = mean(horzcat(sim_B_process_1s,sim_B_process_2s,sim_B_process_3s),2);
sim_C_process = mean(horzcat(sim_C_process_1s,sim_C_process_2s,sim_C_process_3s),2);
F_B_process = (sim_sov_thick.*sim_B_process)*100*1e2*1e-3*0.6; %N
F_C_process = (sim_sov_thick.*sim_C_process)*100*1e2*1e-3*0.6; %N
S=[pi_p ADP_p Pcr_p H_p ATP_p sim_Ftotal_2];
%% plot force
figure(21);clf;
plot(cycles,F_B_process,'linewidth',2,'Color','k'); hold on;
plot(cycles,abs(F_C_process),'linewidth',2,'Color','r'); hold on;box on;
ylim([-50 800]);
xlabel('Cycle index');
ylabel('Force (N)');
set(gca,'Unit','Inches')
p = get(gca,'Position');
set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25]);
exportgraphics(figure(21),fullfile('figure_5_subplots','figure_5_c.pdf'),'BackgroundColor','w','Resolution',300,'ContentType','vector');