%% read the parameters
param_table=readtable('params/params.xlsx');
params = param_table.Values;
m_factor = 0.75:0.025:2;
n = length(m_factor);
%% Simulate force and met profiles for different Pi levels
T_cycles = 60;
force=zeros(T_cycles,n);
Pi_lvl=zeros(T_cycles,n);
ADP_lvl=zeros(T_cycles,n);
PCr_lvl=zeros(T_cycles,n);
H_lvl=zeros(T_cycles,n);
ATP_lvl=zeros(T_cycles,n);
for i=1:n
    params1 = params;
    params1(19) = params(19)*m_factor(i);
    [~,S,~,~,~,~,~]=sim_dynamics(params1);
    force(:,i)=S(:,6);
    Pi_lvl(:,i)=S(:,1);
    ADP_lvl(:,i)=S(:,2);
    PCr_lvl(:,i)=S(:,3);
    H_lvl(:,i)=S(:,4);
    ATP_lvl(:,i)=S(:,5);
end
%% construct figure 5d
figure(62);clf;
plot(m_factor*params(19),force(T_cycles,:),'linewidth',2,'Color','k'); hold on;box on;
xlabel('K_{Pi} (mM)','Interpreter','tex');
ylabel('Force (N)');
ylim([275 315]);
%yticks(280:15:325);
xlim([25 85]);
xticks(25:10:85);
set(gca,'Unit','Inches')
p = get(gca,'Position');
set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25]);
exportgraphics(figure(62),fullfile('figure_5_subplots','figure_5_d.pdf'),'BackgroundColor','w','Resolution',300,'ContentType','vector');