%% add utils to path
addpath([pwd '/amofa_utils']);
addpath([pwd '/ulfmm']);
%% Clustering Demo: Geva's Face
npoints=1000;
pp = [0.2 0.2 0.2 0.2];%0.25
mu1 = [0 0];
mu2 = [0 0];
mu3 = [-1.5 1.5];
mu4 = [1.5 1.5];
mu5= [0 -2];
mu = [mu1' mu2' mu3' mu4' mu5']; % 
covar(:,:,1) = [.01 0; 0 1.25];
covar(:,:,2) = [8 0; 0 8];
covar(:,:,3) = [0.2 0; 0 0.015];
covar(:,:,4) = [0.2 0; 0 0.015];
covar(:,:,5) = [1 0; 0 0.2];
[y,r] = genmixr(npoints,mu,covar,pp);

y=y';
r=r';
[mixture_amofa,training_history]=amofa(y);
figure,
plot_mixture(mixture_amofa',y,r);

plot_evolution3d(training_history,y,r)

[clustering,responsibility,maxlik,params]=  getClustering(y,mixture_amofa);
%% Clustering Demo: Ueda's Spiral
load ueda_spiral
[mixture_amofa_ue,training_history_ue]=amofa(spiral);
plot_mixture3d(mixture_amofa_ue',spiral);

%% Classification Demo: Japanese Phoneme data 
load jpn_phoneme

numC=max(trainlabels);%number of classes

class=cell(numC,1);
priors=zeros(numC,1);
N=size(train,1);
% use parallel processing if available

for i=1:numC
    priors(i)=sum(trainlabels==i)/N;
    [class{i}.mixture, class{i}.trainingHistory_amofa] = amofa(train(trainlabels==i,:));

end


[accuracy_amofa,params_amofa]=scoreMoMoFA(class,test,testlabels,numC);

disp([' Classification accuracy AMoFA ',num2str(accuracy_amofa),' per cent.']);

%% Carry out 10 fold CV test on Japanese Phoneme on data
% NOTE: may take some time
load jpn_phoneme
data=[train;test];
datalabels=[trainlabels;testlabels];
[folds,Performance_Params]=nfoldCVtest(data,datalabels,1,'jpn_phoneme',10);

%%
% Carry out simulations on 100 datasets sampled from 4 Overlapping
% Gaussians (example 2 in the paper)
% NOTE: may take some time
clear all
load 4ovGauss
simulCnt=100;
verbose=0;
covoption=0;
BestKset_AMoFA_4oG=zeros(simulCnt,3);
normInfoDist_4oG=zeros(simulCnt,3);

bestParams_4oG=zeros(simulCnt,3);
cluster_accur_4oG=zeros(1,3);
models=cell(simulCnt,3);
%permute the original dataset to provide a shuffled val set from imofa
perm_set_imo=randperm(1000);
tr_end_imo=700;% for IMOFA training set to val set inflection point
for i=1:simulCnt
    %run AMoFA
    [models{i,1}.mixture,models{i,1}.trainHist_amo] = amofa(squeeze(Y(:,:,i)));
    BestKset_AMoFA_4oG(i,1)=numel(models{i,1}.mixture);

    [models{i,1}.clustering,models{i,1}.responsibility,...
        models{i,1}.maxlik,models{i,1}.params]= ...
        getClustering(squeeze(Y(:,:,i)),models{i,1}.mixture);
    normInfoDist_4oG(i,1)= norm_info_dist(models{i,1}.clustering,squeeze(R(:,i)));
   
    bestParams_4oG(i,1)=models{i,1}.params;
    BestKset_AMoFA_4oG(i,1)=numel(models{i,1}.mixture);
    
    %run IMoFA
    train_imo=squeeze(Y(perm_set_imo(1:tr_end_imo),:,i));
    val_imo=squeeze(Y(perm_set_imo(tr_end_imo+1:end),:,i));
    
    [models{i,2}.mixture,models{i,2}.trainHist_amo] = imofa(train_imo,val_imo);
    
    [models{i,2}.clustering,models{i,2}.responsibility,...
    models{i,2}.maxlik,models{i,2}.params]= ...
    getClustering(squeeze(Y(:,:,i)),models{i,2}.mixture);
    normInfoDist_4oG(i,2)= norm_info_dist(models{i,2}.clustering,squeeze(R(:,i)));
    
    bestParams_4oG(i,2)=models{i,2}.params;
    BestKset_AMoFA_4oG(i,2)=numel(models{i,2}.mixture);
    
     % This part needs ULFMM code, comment out if you did not downlaod it yet
    [models{i,3}.bestk,models{i,3}.bestpp,models{i,3}.bestmu,models{i,3}.bestcov,...
        models{i,3}.dl,models{i,3}.countf] = mixtures4(squeeze(Y(:,:,i))',1,20,0,1e-5,covoption,verbose);
    [models{i,3}.clustering,models{i,3}.responsibility,...
    models{i,3}.maxlik,models{i,3}.params]=  getClusteringGmm(squeeze(Y(:,:,i)),models{i,3}.bestmu,models{i,3}.bestcov,models{i,3}.bestpp);

    BestKset_AMoFA_4oG(i,3)=models{i,3}.bestk;
    normInfoDist_4oG(i,3)= norm_info_dist(models{i,3}.clustering,squeeze(R(:,i)));
    bestParams_4oG(i,3)=models{i,3}.params;

end

mean_NID_4oG=mean(normInfoDist_4oG);
for i=1:3
    cluster_accur_4oG(i)=sum(BestKset_AMoFA_4oG(:,i)==4);
end
str='';
figure,
for i=1:3
   
subplot(1,3,i)
hist(BestKset_AMoFA_4oG(:,i),8)
if (i==1)
    str='AMoFA';
elseif (i==2)
    str='IMoFA';
else
    str='ULFMM';
end
title(str);
axis([1 8 1 simulCnt]) 
end


