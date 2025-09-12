%% pre-simulate iemg data
iemg_table=readtable("../raw_data/iEMG_data_burnley.xlsx");
iemg_val = iemg_table{:,2};
%% reading the parameters 
param_table_dke=readtable('params/params.xlsx');
params_dke = param_table_dke.Values;
%% set n
n = 100;
%% get the range of H+ conc to be evaluated
data_H = readtable('../raw_data/proton_data_broxterman.xlsx'); %H+
H_min = min(data_H{:,2});
H_max = max(data_H{:,2});
del_H = (H_max-H_min)/n;
%% set the range of H+ conc to be evaluated
data_Pi = readtable('../raw_data/Pi_data_broxterman.csv'); %Pi
Pi_min = min(data_Pi{:,2});
Pi_max = max(data_Pi{:,2});
del_Pi = (Pi_max-Pi_min)/n; 
%% create the mesh 
[Pi_mesh, H_mesh] = meshgrid(Pi_min:del_Pi:Pi_max, H_min:del_H:H_max);
%% set the rate of change of all metabolite conc to zero 
met_dyn_set = [0 0 0 0];
%% simulate the force for the diffent H and Pi concentrations in the mesh
S = cell(n+1,n+1);
F_final = zeros(n+1,n+1);
for i=1:n+1
    for j=1:n+1
        Pi_set = Pi_mesh(i,j);
        H_set = H_mesh(i,j);
        [Y,F,~,~,~]=cycle_1_pi_H(params_dke,iemg_val(1),Pi_set,H_set,met_dyn_set);
        S1 = [Y F];
        S{i,j} = S1;
        F_final(i,j)=mean(F);
        clear F S1
    end
end
%% prepare the plot
figure(55)
surf(Pi_mesh, H_mesh*10^6, F_final);
% Add labels and title
xlabel('Pi (mM)');
ylabel('H+ (nM)');
zlabel('Force (N)');

%viewing angle
view(-48,68);

% Add color shading

colormap jet;      % Use 'jet' colormap
colorbar;          % Display colorbar indicating Z values

% Improve appearance
shading interp;    % Smooth shading
% export graph
set(gca,'Unit','Inches')
p = get(gca,'Position');
set(gca,'Unit','Inches','Position',[p(1) p(2) 3 2.25]);
exportgraphics(figure(55),fullfile('figure_4_subplots','figure_4_c.pdf'),'BackgroundColor','w','Resolution',300,'ContentType','vector');