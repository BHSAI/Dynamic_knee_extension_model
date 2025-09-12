%% pre-simulate iemg data
iemg_table=readtable("../raw_data/iEMG_data_burnley.xlsx");
iemg_val = iemg_table{:,2};
%% create the hypothetical parameter sets
param_table_dke=readtable('params/params.xlsx');
params_dke = param_table_dke.Values;
%% set the range of H+ conc to be evaluated
data_Pi = readtable('../raw_data/Pi_data_broxterman.csv'); %Pi
data_resting=readtable("../raw_data/initial_state.xlsx");
Pi_set1 = data_Pi{1,2}:5:data_Pi{60,2};
n=length(Pi_set1);
met_dyn_set = [0 0 0 0];
%% Simulate force and met profiles for different H+ levels
S = cell(n,1);
for i=1:n
    clear P_state;
    [Y,F,~,~,~]=cycle_1_pi(params_dke,iemg_val(1),Pi_set1(i),met_dyn_set);
    P_state = 1-Y(:,1)-Y(:,4)-Y(:,7)-Y(:,10);
    S1 = [Y P_state F];
    S{i} = S1;
end
[w,r] =size(S{1,1});
final_lvl=zeros(n,r);
for i=1:n
    clear S_temp
    S_temp = S{i,1};
    final_lvl(i,:)=S_temp(w,:);
end
%% Plot the change in fractions of states with respect to Pi concentration
figure(32);clf;
y_labels={'N','P','A1','A2','A3'};
column_index = [10 17 1 4 7];
q=length(y_labels);
for i=1:q
    %plot(Pi_set1,final_lvl(:,column_index(i))-final_lvl(1,column_index(i)),'linewidth',2,'DisplayName',y_labels{i}); hold on;
    scatter(Pi_set1,final_lvl(:,column_index(i))-final_lvl(1,column_index(i)),'filled','s','DisplayName',y_labels{i}); hold on; box on;
end
legend('Location','eastoutside','NumColumns',1,'Box','off')
xlim([0 35]);
ylim([-0.1 0.15]);
xticks(0:5:35)
yticks([-0.1 -0.05 0 0.05 0.1 0.15])
ytickformat('%.2f');
xlabel('Pi (mM)');
ylabel('Change in state fractions');
set(gca,'Unit','Inches')
p = get(gca,'Position');
set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25]);
exportgraphics(figure(32),fullfile('figure_4_subplots','figure_4_d.pdf'),'BackgroundColor','w','Resolution',300,'ContentType','vector');