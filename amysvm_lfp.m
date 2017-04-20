% clear
function  [X,Y,ranidx,SVMModel,score1,score2,pfm]=amysvm_lfp(destination)
% destination='i:\160721';
load([destination '\Blksortedcell']);
load([destination '\Amy_unit_block_conLFP.mat']);
load([destination '\Amy_unit_block_oriLFP.mat']);
 %%
%spike based svm 做不出来
%  load([str '\BlkCSF_unsorted']);
% load([str '\BlkNEUF_unsorted']);
% nosorting=1;
% Nsession=2;
% y=repmat([ones(30,1);zeros(90,1)],Nsession,1);
% x=[];
% X=[];
% % figure
% count=0;
% for nblock=1:Nsession:2%numel(Blksortedcell)
%     if sum(sum(Blksortedcell{nblock}))<25;
%         ss=zeros(3,8);
%         for iii=1:numel(Blksortedcell)
%             ss=ss+Blksortedcell{iii};
%         end
%         [b,cellidx]=sort(reshape(ss,24,1),'descend');
%         cellused=zeros(3*8,1);
%         cellused(cellidx(1:12))=1;
%         cellused=reshape(cellused,3,8);
%     else
%         cellused=Blksortedcell{nblock};
%     end
%     for j=1:3
%         for k=1:8
%             if cellused(j,k)
%                 TEMP=[];
%                 for nb=1:Nsession
%                     temp=[];
%                     for nwin=1:size(win,1)
%                         temp=[temp cat(1,mean(BlkCSF{nblock+nb-1}{j,k}(:,win(nwin,1):win(nwin,2)),2),mean(BlkNEUF{nblock+nb-1}{j,k}(:,win(nwin,1):win(nwin,2)),2))];                   
%                     end
%                     TEMP=[TEMP; temp];
%                 end
%                 x=[x TEMP];
% %                 for ntemp=1:size(temp,2)
% %                     hold on
% %                     count=count+1;
% %                     plot(count*5,temp(1:30,ntemp),'ro');
% %                     plot(count*5,mean(temp(1:30,ntemp)),'rx');
% %                     plot(count*5+1,temp(31:120,ntemp),'go');
% %                     plot(count*5+1,mean(temp(31:120,ntemp)),'gx');
% %                 end
%             end
%         end
%     end
%     [coeff,score,latent,tsquared,explained,mu]=pca(x);
%     X=[X;x];
% end
% 
% rearridx=randperm(size(y,1));
% trainidx=1:120;
% testidx=121:240;
% SVMModel=fitcsvm(x(trainidx,:),y(trainidx,:),'KernelScale','auto','Standardize',true,...
%     'OutlierFraction',0.05,'KernelFunction','RBF');
% CVSVMModel = crossval(SVMModel);
% classLoss = kfoldLoss(CVSVMModel);
% [ylabel,scorePred] = kfoldPredict(CVSVMModel);
% label=predict(SVMModel,x(testidx,:));
% sum(label)
% 
% SVMModel=fitcsvm(x,y,'KernelScale','auto','Standardize',true,...
%     'OutlierFraction',0.05,'KernelFunction','RBF');
% CVSVMModel = crossval(SVMModel);
% classLoss = kfoldLoss(CVSVMModel);
% [ylabel,scorePred] = kfoldPredict(CVSVMModel);
% label=predict(SVMModel,x);
% sum(label)

%%
params.Fs=1000;
params.fpass=[10 45];
params.pad=0;
params.tapers=[2 3];
params.trialave=0;
params.err=[1 0.05];

spec=struct;
win1=[401 800];
alllfp=cellfun(@pickuplfp,amyunitblockconLFP,'UniformOutput',false);
[spec.S1,spec.f1]=cellfun(@(A) cellspectrum(A,win1,params),alllfp,'UniformOutput',false);
win2=[801 1200];
params.fpass=[55 90];
[spec.S2,spec.f2]=cellfun(@(A) cellspectrum(A,win2,params),alllfp,'UniformOutput',false);
%%
%cell based,每个cell在每天所有block的lfp做特征分类
% aaa=cell2mat(Blksortedcell);
% bbb=reshape(aaa,3,8,numel(Blksortedcell));
% ccc=mat2cell(bbb,repmat(1,1,3),repmat(1,1,8),numel(Blksortedcell));
% refromedsortedcell=cellfun(@squeeze,ccc,'UniformOutput',false);
% X=cellfun(@unitmergefeature,refromedsortedcell,spec.S1,spec.S2,'UniformOutput',false);
% l=cellfun(@length,X);
% X(l==0)=[];
% refromedsortedcell(l==0)=[];
% Y=repmat([ones(30,1);zeros(90,1)],sum(refromedsortedcell{6}),1);
% id=3;
% ranidx=randperm(size(X{id},2))<floor(size(X{id},2)*0.5)+1;
% SVMModel=fitcsvm(X{id}(:,ranidx)',Y(ranidx),'KernelScale','auto','Standardize',true,...
%     'OutlierFraction',0.1,'KernelFunction','RBF');
% labelb=predict(SVMModel,X{id}(:,ranidx)');
% label=predict(SVMModel,X{id}(:,~ranidx)');
% sum(label)
%%
% block based用每个block中所有有spike的cell的lfp平均后做特征
zzz=cellfun(@(x) cell2mat(x'),spec.S1,'UniformOutput',false);
xxx=cell2mat(zzz(:));
ccc=mat2cell(xxx,size(xxx,1),repmat(size(spec.S1{1,1}{1},2),1,numel(Blksortedcell)));
vvv=cellfun(@(x) mat2cell(x,repmat(size(spec.S1{1,1}{1},1),1,24),size(x,2)),ccc,'UniformOutput',false);
reformedS1=cellfun(@(x) reshape(x,3,8),vvv,'UniformOutput',false);

zzz=cellfun(@(x) cell2mat(x'),spec.S2,'UniformOutput',false);
xxx=cell2mat(zzz(:));
ccc=mat2cell(xxx,size(xxx,1),repmat(size(spec.S2{1,1}{1},2),1,numel(Blksortedcell)));
vvv=cellfun(@(x) mat2cell(x,repmat(size(spec.S2{1,1}{1},1),1,24),size(x,2)),ccc,'UniformOutput',false);
reformedS2=cellfun(@(x) reshape(x,3,8),vvv,'UniformOutput',false);

% X=cellfun(@mergefeature,Blksortedcell,reformedS1,reformedS2,'UniformOutput',false);
% Y=[ones(30,1); zeros(90,1)];
% allX=cell2mat(X)';
% allY=repmat(Y,numel(X),1);
% ranidx=randperm(numel(allY))<numel(allY)*0.5+1;
% SVMModel=fitcsvm(allX(ranidx,:),allY(ranidx),'KernelScale','auto','Standardize',true,...
%     'OutlierFraction',0.05);
% label=predict(SVMModel,allX(~ranidx,:));
% pfm=sum(label==allY(~ranidx,:))

blkX=cellfun(@(a,b,c) merge_cell_blk(a,b,c,30),Blksortedcell,reformedS1,reformedS2,'UniformOutput',false);
matX=cellfun(@(x) cell2mat(x'),blkX,'UniformOutput',false);
X=transpose(cell2mat(matX));
Y=repmat([1;0],size(X,1)/2,1);
ranidx=randperm(numel(Y))<numel(Y)*0.5+1;
SVMModel=fitcsvm(X(ranidx,:),Y(ranidx),'KernelScale','auto','Standardize',true,...
    'OutlierFraction',0.05,'KernelFunction','RBF');
[lb,score1]=predict(SVMModel,X(ranidx,:));
[label,score2]=predict(SVMModel,X(~ranidx,:));
pfm=sum(label==Y(~ranidx))
sprintf([num2str(pfm) '/' num2str(numel(Y(~ranidx)))])
%% 画每个block lfp能量图
% spc1=struct('cs',[],'neu',[]);
% spc2=struct('cs',[],'neu',[]);
% for nblock=1:numel(spec.S1{1,1})
%     for j=1:3
%         for k=1:8
%             if Blksortedcell{nblock}(j,k)
%                 spc1.cs=[sp1.cs spec.S1{j,k}{nblock}(:,1:30)];
%                 spc1.neu=[sp1.neu spec.S1{j,k}{nblock}(:,31:end)];
%                 spc2.cs=[sp1.cs spec.S2{j,k}{nblock}(:,1:30)];
%                 spc2.neu=[sp1.neu spec.S2{j,k}{nblock}(:,31:end)];
%             end
%         end
%     end
% end
% f1=spec.f1{1,1}{1};
% f2=spec.f2{1,1}{1};
% figure
% plot(f1,mean(spc1.cs,2),'r')
% hold on
% plot(f1,mean(spc1.neu,2),'g')
% figure
% plot(f2,mean(spc2.cs,2),'r')
% hold on
% plot(f2,mean(spc2.neu,2),'g')
%% cell based 每个细胞lfp能量图
% f1=spec.f1{1,1}{1};
% f2=spec.f2{1,1}{1};
% for n=1:numel(X)
%     idx1=repmat(1:30,1,floor(size(X{n},2)/120))+repelem(0:floor(size(X{n},2)/120)-1,30)*120;
%     idx2=repmat(31:120,1,floor(size(X{n},2)/120))+repelem(0:floor(size(X{n},2)/120)-1,90)*120;
%     figure
%     subplot(1,2,1)
%     lerb1=std(X{n}(1:18,idx1),1,2)/sqrt(numel(idx1));
%     lerb2=std(X{n}(1:18,idx2),1,2)/sqrt(numel(idx2));
%     errorbar(f1,mean(X{n}(1:18,idx1),2),lerb1,'color','r')
%     hold on
%     errorbar(f1,mean(X{n}(1:18,idx2),2),lerb2,'color','g')
%     subplot(1,2,2)
%     herb1=std(X{n}(19:36,idx1),1,2)/sqrt(numel(idx1));
%     herb2=std(X{n}(19:36,idx2),1,2)/sqrt(numel(idx2));
%     errorbar(f2,mean(X{n}(19:36,idx1),2),herb1,'color','r')
%     hold on
%     errorbar(f2,mean(X{n}(19:36,idx2),2),herb2,'color','g')
% end
function celllfp=pickuplfp(CA)
celllfp=cellfun(@(x) double(cell2mat(x')),CA,'UniformOutput',false);

function [S,f]=cellspectrum(BA,win,param)
[S,f]=cellfun(@(A) mtspectrumc(A(win(1):win(2),:),param),BA,'UniformOutput',false);

function x=mergefeature(Cidx,C1,C2)
s1=C1(Cidx);
ss1=cell2mat(s1');
sss1=reshape(ss1,size(s1{1},1),size(s1{1},2),numel(s1));
s2=C2(Cidx);
ss2=cell2mat(s2');
sss2=reshape(ss2,size(s2{1},1),size(s2{1},2),numel(s2));
x=[mean(sss1,3); mean(sss2,3)];

function ss=merge_cell_blk(Cidx,C1,C2,Ncstrial)
s1=C1(Cidx);
s2=C2(Cidx);
ss=cellfun(@(x,y) [mean(x(:,1:Ncstrial),2) mean(x(:,1+Ncstrial:end),2); mean(y(:,1:Ncstrial),2) mean(y(:,1+Ncstrial:end),2)],...
    s1,s2,'UniformOutput',false);

function x=unitmergefeature(Uidx,C1,C2)
s1=C1(Uidx);
ss1=cell2mat(s1');
s2=C2(Uidx);
ss2=cell2mat(s2');
x=[ss1;ss2];