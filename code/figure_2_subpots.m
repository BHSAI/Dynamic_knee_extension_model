%% create the hypothetical parameter sets
param_table=readtable('params/params.xlsx');
params= param_table.Values;
%experiment results
data_Pcr  = readtable('../raw_data/PCr_data_broxterman.csv'); % Pcr
data_Pi  = readtable('../raw_data/Pi_data_broxterman.csv'); % Pi 
data_ADP  = readtable('../raw_data/ADP_data_broxterman.xlsx'); % ADP 
data_H  = readtable('../raw_data/proton_data_broxterman.xlsx'); % ATP
data_force = readtable('../raw_data/force_data_broxterman.csv'); % force
exp_data = [data_Pi{2:end,2} data_ADP{2:end,2}*10^-3 data_Pcr{2:end,2} data_H{2:end,2} zeros(60,1) data_force{:,2}];
% get simulation results
column_names1={'Pi (mM)','ADP (uM)','PCr (mM)','pH','ATP (mM)','Force (N)'};
column_names2=cell(1,60);
for i=1:60
    column_names2{i}=sprintf('cycle_%d',i);
end
[~,S,sim_Ftotal_cycles,~,~,~,~]=sim_dynamics(params);
dyn_table1=array2table(S,"VariableNames",column_names1);
dyn_table2=array2table(sim_Ftotal_cycles,"VariableNames",column_names2);
% write tables
writetable(dyn_table1,fullfile('figure_2_subplots/','dynamics1.xlsx'));
writetable(dyn_table2,fullfile('figure_2_subplots/','dynamics2.xlsx'));
cycles = 1:60;
%plot simulations and experiment results
labels_y={'Pi (mM)','ADP (uM)','PCr (mM)','pH','ATP (mM)','Force (N)'};
filename={'Pi_fit.pdf','ADP_fit.pdf','PCr_fit.pdf','H_fit.pdf','ATP_fit.pdf','force_fit.pdf'};
[~,q]=size(S);
rmsd = zeros (q,1);
rmsd_unit = {'mM','uM','mM','pH','mM','N'};
for i=1:q
    figure(50+i);clf;
    x_txt=0.99;
    y_txt=0.05;
    if i==4
        plot(cycles,-log10(S(:,i)*10^-3),'linewidth',2,'Color','k'); hold on;
        scatter(cycles,-log10(exp_data(:,i)*10^-3),20,'MarkerEdgeColor',[1 1 1]*0.5,'MarkerFaceColor',[1 1 1]*0.5,'Marker','square');
        %diff = (S(:,i) - exp_data(:,i))*10^6;
        diff = (-log10(S(:,i)*10^-3))-(-log10(exp_data(:,i)*10^-3));
        rmsd(i)= rms(diff);
        rmse_string = sprintf('RMSE: %1.1f',rmsd(i));
        text(x_txt,y_txt,rmse_string,"Units","normalized","HorizontalAlignment","right","VerticalAlignment","baseline","FontSize",8);
        ylim([6.0 7.5]);
        ytickformat('%.1f')
    elseif i==5
        plot(cycles,S(:,i),'linewidth',2,'Color','k');ylim([0 inf]);
    else
        plot(cycles,S(:,i),'linewidth',2,'Color','k'); hold on;
        scatter(cycles,exp_data(:,i),20,'MarkerEdgeColor',[1 1 1]*0.5,'MarkerFaceColor',[1 1 1]*0.5,'Marker','square');
        diff = S(:,i) - exp_data(:,i);
        rmsd(i)= rms(diff);
        rmse_string = sprintf('RMSE: %1.1f %s',rmsd(i),rmsd_unit{i});
        text(x_txt,y_txt,rmse_string,"Units","normalized","HorizontalAlignment","right","VerticalAlignment","baseline","FontSize",8);
        if i==1
            ylim([0 40]);
        elseif i==3
            ylim([0 25]);
        elseif i==6
            ylim([0 inf])
        end
    end
    ylabel(labels_y(i),'fontsize',11)
    xlabel('Cycle index','fontsize',11)
    set(gca,'Unit','Inches')
    p = get(gca,'Position');
    set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25]);
    exportgraphics(figure(50+i),fullfile('figure_2_subplots/',filename{i}),'BackgroundColor','w','Resolution',300,'ContentType','vector');
end
%compile the rmsd
t_rmsd=table(rmsd,'RowNames',labels_y');
writetable(t_rmsd,fullfile('figure_2_subplots/','rmsd.xlsx'),'WriteRowNames',true);