%% Figure2------Heatmap
clear all

%% load data 
load('HCP_280_IQ.mat')
load('spearman.mat')
IQ_titles = insertBefore(IQ_titles,'_','\');



%% characteristic path length + edge density

binL = zeros(size(netsAll,3),1);
weiL = zeros(size(netsAll,3),1);
metL = zeros(size(netsAll,3),1);
density = zeros(size(netsAll,3),1);
for n = 1:  size(netsAll,3)
    %binary
    binD = distance_bin(netsAll(:,:,n) ~=0);
    binL(n) = charpath(binD);
    
    %weighted  ---the inverse of the number of streamlines is used
    weiD = distance_wei(1./netsAll(:,:,n));
    weiL(n) = charpath(weiD);
    
    %metric   ---use e_distsAll
    metD = distance_wei(e_distsAll(:,:,n).* (netsAll(:,:,n)~=0));
    netsAll(:,:,n)~=0
    metL(n) = charpath(metD);
    
    %edge density
    density(n) = density_und(netsAll(:,:,n));
     
end




%% clustering coefficient:
binC = zeros(size(netsAll,3),1);
weiC = zeros(size(netsAll,3),1);

for n = 1: size(netsAll,3)
    %binary
    binC(n) = mean(clustering_coef_bu(netsAll(:,:,n) ~=0));
    
    %weighted
    net_nrm = weight_conversion(netsAll(:,:,n), 'normalize');
    weiC(n) = mean(clustering_coef_wu(net_nrm));  
end


%% e_dist_bc 
EDISTbc_all = zeros(size(netsAll,3),1);
EDISTbc_long = zeros(size(netsAll,3),1);
EDISTbc_short = zeros(size(netsAll,3),1);


for n = 1: size(netsAll,3)
    tic;
    tmp_net = netsAll(:,:,n);
    tmp_e_dist = e_distsAll(:,:,n);
    bin_net = (tmp_net~= 0);

    e_dist_edge = tmp_e_dist.*(triu(tmp_net)~= 0);
    tmp = e_dist_edge(e_dist_edge ~=0);
    m_edge = median(tmp);
    
    e_dist_edge = tmp_e_dist.*(tmp_net~= 0);
    e_dist_edge_long = (tmp_e_dist >= m_edge).*(tmp_net~= 0);
    e_dist_edge_short = (tmp_e_dist < m_edge).*(tmp_net~= 0);
    
    
    edistbc_all2 = 0;
    
    pair_count_all2 = 0;
    
    edistbc_long2 = 0;
    
    pair_count_long2 = 0;
    
    edistbc_short2 = 0;
    
    pair_count_short2 = 0;
    
    for i = 1: size(tmp_net,1)
        connect_nodes = find(e_dist_edge(i,:)~=0);
        for j = 1: length(connect_nodes)
            for k = j+1 : length(connect_nodes)
                a = tmp_e_dist(connect_nodes(j),connect_nodes(k));
                b = bin_net(connect_nodes(j),connect_nodes(k));
                
                edistbc_all2 = edistbc_all2 + a.* (1 - b);
                
                pair_count_all2 = pair_count_all2 + 1*(1 - b);
            end
        end
        
        connect_nodes = find(e_dist_edge_long(i,:)~=0);
        for j = 1: length(connect_nodes)
            for k = j+1 : length(connect_nodes)
                a = tmp_e_dist(connect_nodes(j),connect_nodes(k));
                b = bin_net(connect_nodes(j),connect_nodes(k));
                
                edistbc_long2 = edistbc_long2 + a.* (1 - b);
                
                pair_count_long2 = pair_count_long2 + 1*(1 - b);
            end
        end
        
        connect_nodes = find(e_dist_edge_short(i,:)~=0);
        for j = 1: length(connect_nodes)
            for k = j+1 : length(connect_nodes)
                a = tmp_e_dist(connect_nodes(j),connect_nodes(k));
                b = bin_net(connect_nodes(j),connect_nodes(k));
                
                edistbc_short2 = edistbc_short2 + a.* (1 - b);
                
                pair_count_short2 = pair_count_short2 + 1*(1 - b);
            end
        end
    end
    
    EDISTbc_all(n) = edistbc_all2 / pair_count_all2;
    EDISTbc_long(n) = edistbc_long2 /pair_count_long2;
    EDISTbc_short(n) = edistbc_short2 /pair_count_short2;
    toc;
    
end


%% Volume convex hull and Rentian Exponent
V_convex_hull = zeros(size(netsAll,3),1);
for n = 1:  size(netsAll,3)
    [~,V_convex_hull(n)] = convhull(xyzAll(:,1,n),xyzAll(:,2,n),xyzAll(:,3,n));


end



%% save data
save('statistic.mat', 'V_convex_hull', 'EDISTbc_all', 'EDISTbc_long', 'EDISTbc_short', 'metL', 'binC', 'weiC', 'density', 'binL', 'weiL' )


%% integrate data
spatial_features = [V_convex_hull EDISTbc_all EDISTbc_long EDISTbc_short metL];
spatial_features_titles = {'Convex Hull Volume' 'Dist.BC (global)' 'Dist.BC (long)' 'Dist.BC (short)' 'L_m'};

topo_features = [binC weiC density binL weiL];
topo_features_titles = {'C_u' 'C_w' 'Edge Density' 'L_u' 'L_w'};

class2_titles = IQ_titles([49, 12, 10, 3, 5, 14, 30, 7, 1, 43, 41]);
class2_IQ = IQ(:,[49, 12, 10, 3, 5, 14, 30, 7, 1, 43, 41]);

frequency_titles={'40 Hz','20 Hz','10 Hz'};

figure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[RHO_topo_IQ2,PVAL_topo_IQ2] = corr(topo_features,class2_IQ ,'Type','Spearman');

subplot(3,5,1);
h_topo2_0 = heatmap(class2_titles(1), topo_features_titles, RHO_topo_IQ2(:,1), 'Colormap',redblue(64),'ColorbarVisible', 'off', 'GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_topo2_0.XDisplayLabels = nan(size(h_topo2_0.XDisplayData));
h_topo2_0.Title='Full IQ';
set(gca,'FontSize',9);

subplot(3,5,2);
h_topo2_1 = heatmap(class2_titles(2:3), topo_features_titles, RHO_topo_IQ2(:,2:3), 'Colormap',redblue(64), 'ColorbarVisible', 'off','GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_topo2_1.XDisplayLabels = nan(size(h_topo2_1.XDisplayData));
h_topo2_1.YDisplayLabels = nan(size(h_topo2_1.YDisplayData));
h_topo2_1.Title='Crystallized';
set(gca,'FontSize',9);

subplot(3,5,3);
h_topo2_2 = heatmap(class2_titles(4:6), topo_features_titles, RHO_topo_IQ2(:,4:6), 'Colormap',redblue(64), 'ColorbarVisible', 'off','GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_topo2_2.XDisplayLabels = nan(size(h_topo2_2.XDisplayData));
h_topo2_2.YDisplayLabels = nan(size(h_topo2_2.YDisplayData));
h_topo2_2.Title='Processing speed';
set(gca,'FontSize',9);

subplot(3,5,4);
h_topo2_3 = heatmap(class2_titles(7:8), topo_features_titles, RHO_topo_IQ2(:,7:8), 'Colormap',redblue(64), 'ColorbarVisible', 'off','GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_topo2_3.XDisplayLabels = nan(size(h_topo2_3.XDisplayData));
h_topo2_3.YDisplayLabels = nan(size(h_topo2_3.YDisplayData));
h_topo2_3.Title='Visuospatial ability';
set(gca,'FontSize',9);

subplot(3,5,5);
h_topo2_4 = heatmap(class2_titles(9:11), topo_features_titles, RHO_topo_IQ2(:,9:11), 'Colormap',redblue(64), 'ColorbarVisible', 'on','GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_topo2_4.XDisplayLabels = nan(size(h_topo2_4.XDisplayData));
h_topo2_4.YDisplayLabels = nan(size(h_topo2_4.YDisplayData));
h_topo2_4.Title='Memory';
set(gca,'FontSize',9);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[RHO_spatial_IQ2,PVAL_spatial_IQ2] = corr(spatial_features,class2_IQ ,'Type','Spearman');

subplot(3,5,6); 
h_spatial2_0 = heatmap(class2_titles(1), spatial_features_titles, RHO_spatial_IQ2(:,1), 'Colormap',redblue(64), 'ColorbarVisible', 'off','GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_spatial2_0.XDisplayLabels = nan(size(h_spatial2_0.XDisplayData));
h_spatial2_0.Title='Full IQ';
set(gca,'FontSize',9);


subplot(3,5,7);
h_spatial2_1 = heatmap(class2_titles(2:3), spatial_features_titles, RHO_spatial_IQ2(:,2:3), 'Colormap',redblue(64), 'ColorbarVisible', 'off','GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_spatial2_1.XDisplayLabels = nan(size(h_spatial2_1.XDisplayData));
h_spatial2_1.YDisplayLabels = nan(size(h_spatial2_1.YDisplayData));
h_spatial2_1.Title='Crystallized';
set(gca,'FontSize',9);


subplot(3,5,8);
h_spatial2_2 = heatmap(class2_titles(4:6), spatial_features_titles, RHO_spatial_IQ2(:,4:6), 'Colormap',redblue(64), 'ColorbarVisible', 'off','GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_spatial2_2.XDisplayLabels = nan(size(h_spatial2_2.XDisplayData));
h_spatial2_2.YDisplayLabels = nan(size(h_spatial2_2.YDisplayData));
h_spatial2_2.Title='Processing speed';
set(gca,'FontSize',9);


subplot(3,5,9);
h_spatial2_3 = heatmap(class2_titles(7:8), spatial_features_titles, RHO_spatial_IQ2(:,7:8), 'Colormap',redblue(64), 'ColorbarVisible', 'off','GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_spatial2_3.XDisplayLabels = nan(size(h_spatial2_3.XDisplayData));
h_spatial2_3.YDisplayLabels = nan(size(h_spatial2_3.YDisplayData));
h_spatial2_3.Title='Visuospatial ability';
set(gca,'FontSize',9);


subplot(3,5,10);
h_spatial2_4 = heatmap(class2_titles(9:11), spatial_features_titles, RHO_spatial_IQ2(:,9:11), 'Colormap',redblue(64),'ColorbarVisible', 'on', 'GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_spatial2_4.XDisplayLabels = nan(size(h_spatial2_4.XDisplayData));
h_spatial2_4.YDisplayLabels = nan(size(h_spatial2_4.YDisplayData));
h_spatial2_4.Title='Memory';
set(gca,'FontSize',9);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,5,11);
h_frequency3_0 = heatmap(class2_titles(1), frequency_titles, spearman(:,1), 'Colormap',redblue(64), 'ColorbarVisible', 'off', 'GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_frequency3_0.Title='Full IQ';
set(gca,'FontSize',9);


subplot(3,5,12);
h_frequency3_1 = heatmap(class2_titles(2:3), frequency_titles, spearman(:,2:3), 'Colormap',redblue(64), 'ColorbarVisible', 'off','GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_frequency3_1.YDisplayLabels = nan(size(h_frequency3_1.YDisplayData));
h_frequency3_1.Title='Crystallized';
set(gca,'FontSize',9);


subplot(3,5,13);
h_frequency3_2 = heatmap(class2_titles(4:6), frequency_titles, spearman(:,4:6), 'Colormap',redblue(64), 'ColorbarVisible', 'off','GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_frequency3_2.YDisplayLabels = nan(size(h_frequency3_2.YDisplayData));
h_frequency3_2.Title='Processing speed';
set(gca,'FontSize',9);


subplot(3,5,14);
h_frequency3_3 = heatmap(class2_titles(7:8), frequency_titles, spearman(:,7:8), 'Colormap',redblue(64), 'ColorbarVisible', 'off','GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_frequency3_3.YDisplayLabels = nan(size(h_frequency3_3.YDisplayData));
h_frequency3_3.Title='Visuospatial ability';
set(gca,'FontSize',9);


subplot(3,5,15);
h_frequency3_4 = heatmap(class2_titles(9:11), frequency_titles, spearman(:,9:11), 'Colormap',redblue(64), 'GridVisible', 'on', 'CellLabelColor','none', 'fontsize',11);
caxis([-0.3, 0.3]);
h_frequency3_4.YDisplayLabels = nan(size(h_frequency3_4.YDisplayData));
h_frequency3_4.Title='Memory';
set(gca,'FontSize',9);


set(gcf,'Color','w');




