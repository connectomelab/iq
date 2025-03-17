%% Figure S2------to obtain the optimal pair of K and V
%The first 280 subjects are all and the following 15 subjects are just an example

%% 280 subjects
clc;clear all
load('R_SmeanSync4015.mat');load('R_SmetaStab4015.mat');load('MSE4015end.mat');
load('R_SmeanSync4030.mat');load('R_SmetaStab4030.mat');load('MSE4030end.mat');
load('R_SmeanSync4045.mat');load('R_SmetaStab4045.mat');load('MSE4045end.mat');
load('R_SmeanSync4060.mat');load('R_SmetaStab4060.mat');load('MSE4060end.mat');
load('R_SmeanSync4075.mat');load('R_SmetaStab4075.mat');load('MSE4075end.mat');
load('R_SmeanSync4090.mat');load('R_SmetaStab4090.mat');load('MSE4090end.mat');
load('R_SmeanSync40105.mat');load('R_SmetaStab40105.mat');load('MSE40105end.mat');
load('R_SmeanSync40120.mat');load('R_SmetaStab40120.mat');load('MSE40120end.mat');
load('R_SmeanSync40135.mat');load('R_SmetaStab40135.mat');load('MSE40135end.mat');
load('R_SmeanSync40150.mat');load('R_SmetaStab40150.mat');load('MSE40150end.mat');
load('R_SmeanSync40165.mat');load('R_SmetaStab40165.mat');load('MSE40165end.mat');
load('R_SmeanSync40180.mat');load('R_SmetaStab40180.mat');load('MSE40180end.mat');
load('R_SmeanSync40195.mat');load('R_SmetaStab40195.mat');load('MSE40195end.mat');
load('R_SmeanSync40210.mat');load('R_SmetaStab40210.mat');load('MSE40210end.mat');
load('R_SmeanSync40225.mat');load('R_SmetaStab40225.mat');load('MSE40225end.mat');
load('R_SmeanSync40240.mat');load('R_SmetaStab40240.mat');load('MSE40240end.mat');
load('R_SmeanSync40255.mat');load('R_SmetaStab40255.mat');load('MSE40255end.mat');
load('R_SmeanSync40270.mat');load('R_SmetaStab40270.mat');load('MSE40270end.mat');
load('R_SmeanSync40280.mat');load('R_SmetaStab40280.mat');load('MSE40280end.mat');


SumR_SmeanSync40 = R_SmeanSync4015+R_SmeanSync4030+R_SmeanSync4045+R_SmeanSync4060+R_SmeanSync4075+R_SmeanSync4090+R_SmeanSync40105+R_SmeanSync40120+R_SmeanSync40135+R_SmeanSync40150+R_SmeanSync40165+R_SmeanSync40180+R_SmeanSync40195+R_SmeanSync40210+R_SmeanSync40225+R_SmeanSync40240+R_SmeanSync40255+R_SmeanSync40270+R_SmeanSync40280;
SumR_SmetaStab40 = R_SmetaStab4015+R_SmetaStab4030+R_SmetaStab4045+R_SmetaStab4060+R_SmetaStab4075+R_SmetaStab4090+R_SmetaStab40105+R_SmetaStab40120+R_SmetaStab40135+R_SmetaStab40150+R_SmetaStab40165+R_SmetaStab40180+R_SmetaStab40195+R_SmetaStab40210+R_SmetaStab40225+R_SmetaStab40240+R_SmetaStab40255+R_SmetaStab40270+R_SmetaStab40280;
Sum_MSE40end = MSE4015end+MSE4030end+MSE4045end+MSE4060end+MSE4075end+MSE4090end+MSE40105end+MSE40120end+MSE40135end+MSE40150end+MSE40165end+MSE40180end+MSE40195end+MSE40210end+MSE40225end+MSE40240end+MSE40255end+MSE40270end+MSE40280end;

R_SmeanSync40_final = (1/19)*SumR_SmeanSync40;
R_SmetaStab40_final = (1/19)*SumR_SmetaStab40;
MSE40end_final = (1/19)*Sum_MSE40end;
R_SmetaXMSE40 = R_SmetaStab40_final .* (1 - MSE40end_final);



%%min-max normalization
Min_meta = min(min(R_SmetaStab40_final));
Max_meta = max(max(R_SmetaStab40_final));

Min_MSE = min(min(MSE40end_final));
Max_MSE = max(max(MSE40end_final));

R_SmetaStab40_norm = (R_SmetaStab40_final - Min_meta)/(Max_meta - Min_meta);

R_SMSE40_norm = (MSE40end_final - Min_MSE)/(Max_MSE - Min_MSE);

R_SmetaXMSE40_norm = R_SmetaStab40_norm .* (1 - R_SMSE40_norm);

%%imagesc
x1=linspace(0.1,30,30);
y1=linspace(0.1,30,30);
figure(1)
subplot(2,2,1);
imagesc(x1,y1,R_SmeanSync40_final)
xticks([0.1 10 20 30]);
yticks([0.1 10 20 30]);
xlabel('Conduction Velocity');
ylabel('Global Coupling'); 
title('Synchrony');
colorbar;

subplot(2,2,2);
imagesc(x1,y1,R_SmetaStab40_final);
xticks([0.1 10 20 30]);
yticks([0.1 10 20 30]);
xlabel('Conduction Velocity');
ylabel('Global Coupling');
title('Metastability');
colorbar;

% MSE
subplot(2,2,3);
imagesc(x1,y1,MSE40end_final);
xticks([0.1 10 20 30]);
yticks([0.1 10 20 30]);
xlabel('Conduction Velocity');
ylabel('Global Coupling');
title('MSE');
colorbar;

% metastability*(1-MSE)
subplot(2,2,4);
imagesc(x1,y1,R_SmetaXMSE40_norm);
xticks([0.1 10 20 30]);
yticks([0.1 10 20 30]);
xlabel('Conduction Velocity');
ylabel('Global Coupling');
title('metastability*(1-MSE)');
colorbar;


%% 15 subjects

%%min-max normalization
% Min_meta = min(min(R_SmetaStab4015));
% Max_meta = max(max(R_SmetaStab4015));
% Min_MSE = min(min(MSE4015end));
% Max_MSE = max(max(MSE4015end));
% R_SmetaStab4015_norm = (R_SmetaStab4015 - Min_meta)/(Max_meta - Min_meta);
% 
% R_SMSE4015_norm = (MSE4015end - Min_MSE)/(Max_MSE - Min_MSE);
% R_SmetaXMSE4015_norm = R_SmetaStab4015_norm .* (1 - R_SMSE4015_norm);
% 
% %%imagesc
% 
% x1=linspace(0.1,30,30);
% y1=linspace(0.1,30,30);
% 
% figure(1)
% subplot(2,2,1);
% imagesc(x1,y1,R_SmeanSync4015)
% xticks([0.1 10 20 30]);
% yticks([0.1 10 20 30]);
% xlabel('Conduction Speed');
% ylabel('Global Coupling'); 
% title('Synchrony mean(KOP)');
% colorbar;
% 
% subplot(2,2,2);
% imagesc(x1,y1,R_SmetaStab4015);
% xticks([0.1 10 20 30]);
% yticks([0.1 10 20 30]);
% xlabel('Conduction Speed');
% ylabel('Global Coupling');
% title('Metastability std(KOP)');
% colorbar;
% 
% % MSE
% subplot(2,2,3);
% imagesc(x1,y1,MSE4015end);
% xticks([0.1 10 20 30]);
% yticks([0.1 10 20 30]);
% xlabel('Conduction Speed');
% ylabel('Global Coupling');
% title('MSE');
% colorbar;
% 
% % metastability*(1-MSE)
% subplot(2,2,4);
% imagesc(x1,y1,R_SmetaXMSE4015_norm);
% xticks([0.1 10 20 30]);
% yticks([0.1 10 20 30]);
% xlabel('Conduction Speed');
% ylabel('Global Coupling');
% title('metastability*(1-MSE)');
% colorbar;
