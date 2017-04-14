% clear
function [X,SVMModel]=amysvm_spk()
win=[201 600;601 950];
str='i:\160720';
load([str '\Blksortedcell']);
load([str '\BlkCSF_unsorted']);
load([str '\BlkNEUF_unsorted']);
%% 旧版 做不出来
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
%% 新版 向量化 block based
ss=zeros(3,8);
Nused=20;
for iii=1:numel(Blksortedcell)
    ss=ss+Blksortedcell{iii};
end
[b,cellidx]=sort(reshape(ss,24,1),'descend');
cellused=zeros(3*8,1);
cellused(cellidx(1:Nused))=true;
cellused=reshape(cellused,3,8);
blkft=cellfun(@(x,y) block_merge_feature(x,y,win,cellused),BlkCSF,BlkNEUF,'UniformOutput',false);
X=cell2mat(blkft');
Y=repmat([1 1 1 0 0 0],1,size(X,1)/6);
ranidx=randperm(numel(Y))<numel(Y)/2+1;
SVMModel=fitcsvm(X(ranidx,:),Y(ranidx),'KernelScale','auto','Standardize',true,...
    'OutlierFraction',0.05,'KernelFunction','RBF');
label=predict(SVMModel,X(~ranidx,:));
pfm=sum(label==Y(~ranidx)')
sprintf([num2str(pfm) '/' num2str(numel(Y(~ranidx)))])
function [c]=unit_merge_feature(a,b,win)
aft1=arrayfun(@(y,z) mean(a(:,y:z),2),win(:,1),win(:,2),'UniformOutput',false);
aft1a=cell2mat(aft1');
aft1m=[mean(aft1a(1:size(aft1a,1)/3,:));mean(aft1a(1+size(aft1a,1)/3:2*size(aft1a,1)/3,:));mean(aft1a(1+2*size(aft1a,1)/3:end,:))];
aft2=arrayfun(@(y,z) mean(b(:,y:z),2),win(:,1),win(:,2),'UniformOutput',false);
aft2a=cell2mat(aft2');
aft2m=[mean(aft2a(1:size(aft2a,1)/3,:));mean(aft2a(1+size(aft2a,1)/3:2*size(aft2a,1)/3,:));mean(aft2a(1+2*size(aft2a,1)/3:end,:))];
c=[aft1m; aft2m];

function [cc]=block_merge_feature(a,b,win,cellused)
fc=cellfun(@(x,y) unit_merge_feature(x,y,win),a(logical(cellused)),b(logical(cellused)),'UniformOutput',false);
cc=cell2mat(fc');

