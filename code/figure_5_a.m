%% sort sensitivities
params_table = readtable('params/params.xlsx');
params = params_table.Parameter;
sens_force = readtable('../raw_data/Local_sens_Force.xlsx');
[~,I]=sort(abs(sens_force.c60),'descend');
sens_c60_sort=sens_force.c60(I);
%% plot sensitivities
figure(7);clf;
scatter(1:1:23,sens_c60_sort,30,"blue","square","filled"); hold on;box on;
ylim([-2 2]);
xlabel('Parameters');
ylabel('Sensitivity');
set(gca,'Unit','Inches')
yline(0,'-.');
ax=gca;
ax.XTick = 1:1:24;
ax.XTickLabels = params(I);
ax.TickLabelInterpreter = "tex";
set(gca,'Unit','Inches')
p = get(gca,'Position');
set(gca,'Unit','Inches','Position',[p(1)-0.4 p(2) 5.4 1.25]);
exportgraphics(figure(7),fullfile(pwd,'figure_5_subplots','figure_5_a.pdf'),'BackgroundColor','w','Resolution',300,'ContentType','vector');