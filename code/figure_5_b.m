%% sort sensitivities in the same order as local sensitivities
params_table = readtable('params/params.xlsx');
params = params_table.Parameter;
sens_force = readtable('../raw_data/Local_sens_Force.xlsx');
[~,I]=sort(abs(sens_force.c60),'descend');
sens = readmatrix('../raw_data/global_sens_cyc_10000_0.1.xlsx');
sens_sort = sens(I,:);
%% plot sensitivities
figure(8);clf;
boxplot(sens_sort','Colors','b',"Symbol",".","OutlierSize",0.1);hold on;box on;
ylim([-2,2]);
xlabel('Parameters');
ylabel('Sensitivity');
yline(0,'-.');
ax=gca;
ax.XTickLabels = params(I);
ax.TickLabelInterpreter = "tex";
set(gca,'Unit','Inches')
p = get(gca,'Position');
set(gca,'Unit','Inches','Position',[p(1)-0.4 p(2) 5.4 1.25]);
exportgraphics(figure(8),fullfile(pwd,'figure_5_subplots','figure_5_b.pdf'),'BackgroundColor','w','Resolution',300,'ContentType','vector');