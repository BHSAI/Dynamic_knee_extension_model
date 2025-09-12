%% pre-simulate iemg data
iemg_table=readtable("../raw_data/iEMG_data_burnley.xlsx");
iemg_time = iemg_table{:,1}; iemg_val = iemg_table{:,2};
[p_iemg,~,~]=poly2(iemg_time,iemg_val);
cycle_index_exp=1:60; 
cycle_start_time=1:5:300;
cycles=1:1:max(cycle_index_exp); 
m=length(cycles);
T_cycles=(m+(m*0.10));
iemg_profile = zeros(1,T_cycles);
for i=1:m
    iemg_profile(1,i)=((p_iemg(1)*(cycle_start_time(i)^2))+(p_iemg(2)*(cycle_start_time(i)^1))+p_iemg(3))/100;
end
for i=m+1:T_cycles
    iemg_profile(i)=iemg_profile(m);
end
%% set the inputs for the model
param_table=readtable('params/params.xlsx');
params = param_table.Values;
cut_off=m;
Pi_set = 1:0.01:2;
n=length(Pi_set);
met_dyn_set = [0 0 0 0];
%% Simulate force and met profiles for different Pi levels
force=zeros(T_cycles,n);
Pi_lvl=zeros(T_cycles,n);
ADP_lvl=zeros(T_cycles,n);
PCr_lvl=zeros(T_cycles,n);
H_lvl=zeros(T_cycles,n);
ATP_lvl=zeros(T_cycles,n);
for i=1:n
    iemg_temp=iemg_profile;
    [~,S,~,~,~,~,~]=eval_Pi(params,iemg_profile,T_cycles, cut_off, Pi_set(i), met_dyn_set);
    force(:,i)=S(:,6);
    Pi_lvl(:,i)=S(:,1);
    ADP_lvl(:,i)=S(:,2);
    PCr_lvl(:,i)=S(:,3);
    H_lvl(:,i)=S(:,4);
    ATP_lvl(:,i)=S(:,5);
end
%% construct figure 4a
figure(4);clf;
plot(Pi_set,force(T_cycles,:)/force(m,1),'linewidth',2,'Color','k'); hold on;box on;
ylim([0 2]);
ytickformat('%.1f');
xtickformat('%.1f');
xlabel('Pi (Normalized)');
ylabel('Force (Normalized)');
set(gca,'Unit','Inches')
p = get(gca,'Position');
set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25],'YGrid','on');
exportgraphics(figure(4),fullfile('figure_4_subplots','figure_4_a.pdf'),'BackgroundColor','w','Resolution',300,'ContentType','vector');