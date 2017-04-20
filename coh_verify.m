clear
folderlist={
    'I:\160712';
    'I:\160713';
    'I:\160714';
    'I:\160715';
    'I:\160716';
    'I:\160719';
    'I:\160720';
    'I:\160721';
    'I:\160722';
    };
% for n=1:numel(folderlist)
%     n,
% %     MergeV1LFP(folderlist{n},'R',[-0.3 1.1],'ORI');
% %     [LFPCoh1]=V1_Amy_lfpconherency(folderlist{n},[301 1300],1);
%     [LFPCoh2]=V1_Amy_lfpconherency(folderlist{n},[401 800],1);
% %     [LFPCoh2]=V1_Amy_lfpconherency(folderlist{n},[501 800],1);
% %     [LFPCoh3]=V1_Amy_lfpconherency(folderlist{n},[901 1300],1);
% end
% sprintf('done')
%%
%%
stimC1=[];
stimC2=[];
for n=1:numel(folderlist)
    load([folderlist{n} '\LFPCoh101-500ms.mat'])
    stimf=LFPCoh{1}{1}.f1{1};
    for i=1:numel(LFPCoh)
        for j=1:numel(LFPCoh{i})
            stimC1=[stimC1 LFPCoh{i}{j}.C1];
            stimC2=[stimC2 LFPCoh{i}{j}.C2];
        end
    end
end

figure
plot(stimf(1:18),mean(stimC1(1:18,:),2),'r')
hold on
plot(stimf(1:18),mean(stimC2(1:18,:),2),'g')
title('stim time Coh')

alpha1=sum(stimC1(3:3,:),1);
alpha2=sum(stimC2(3:3,:),1);
alphadiff=alpha1-alpha2;
[~,alphaidx]=sort(alphadiff,'ascend');
beta1=sum(stimC1(10:10,:),1);
beta2=sum(stimC2(10:10,:),1);
betadiff=beta1-beta2;
[~,betaidx]=sort(betadiff,'descend');

figure
plot(stimf(1:18),mean(stimC1(1:18,betaidx(1:5000)),2),'r')
hold on
plot(stimf(1:18),mean(stimC2(1:18,betaidx(1:5000)),2),'g')
title('stim time Coh')

figure
plot(stimf(1:18),mean(stimC1(1:18,alphaidx(1:5000)),2),'r')
hold on
plot(stimf(1:18),mean(stimC2(1:18,alphaidx(1:5000)),2),'g')
title('stim time Coh')

figure
plot(alphadiff,betadiff,'r.')
