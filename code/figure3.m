%% simulate iemg data sets
iemg_table=readtable("../raw_data/iEMG_data_burnley.xlsx");
iemg_time = iemg_table{:,1}; iemg_val = iemg_table{:,2};
[p_iemg,~,~]=poly2(iemg_time,iemg_val);
cycle_index_exp=1:60; 
cycle_start_time=1:5:300;
cycles=1:1:max(cycle_index_exp); 
m=length(cycles);
iemg_profile = zeros(2,m);
for i=1:m
    iemg_profile(1,i)=((p_iemg(1)*(cycle_start_time(i)^2))+(p_iemg(2)*(cycle_start_time(i)^1))+p_iemg(3))/100;
    iemg_profile(2,i)=((p_iemg(1)*(cycle_start_time(1)^2))+(p_iemg(2)*(cycle_start_time(1)^1))+p_iemg(3))/100;
end
%% Simulate DKE cycles for the different iemg datasets
clear i
param_table=readtable('params/params.xlsx');
params = param_table.Values;
n= size(iemg_profile,1);
force=zeros(n,m);
for i=1:n
    clear iemg_temp
    iemg_temp=iemg_profile(i,:);
    [~,S,~,~,~,~,~]=eval_emg(params,iemg_temp,m);
    force(i,:)=S(:,6);
end
%% Construct figure 3b
figure(2);clf;
x=1:m;
plot(x,force(1,:),'linewidth',2,'Color','k'); hold on;
plot(x,force(2,:),'linewidth',2,'Color','k','LineStyle','--');
xlim([0 60]);
ylim([200 700]);
xlabel('Cycle index');
ylabel('Force (N)');
set(gca,'Unit','Inches')
p = get(gca,'Position');
set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25]);
exportgraphics(figure(2),fullfile('figure_3_subplots','figure_6_b.pdf'),'BackgroundColor','w','Resolution',300,'ContentType','vector');
%% Construct figure 3a
figure(3);clf;
x=1:m;
plot(x,iemg_profile(1,:)*100,'linewidth',2,'Color','k'); hold on;
plot(x,iemg_profile(2,:)*100,'linewidth',2,'Color','k','LineStyle','--');
xlim([0 60]);
ylim([40 120]);
xlabel('Cycle index');
ylabel('iEMG (%)');
set(gca,'Unit','Inches')
p = get(gca,'Position');
set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25]);
exportgraphics(figure(3),fullfile('figure_3_subplots','figure_6_a.pdf'),'BackgroundColor','w','Resolution',300,'ContentType','vector');