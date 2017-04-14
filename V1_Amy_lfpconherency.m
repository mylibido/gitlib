function [LFPCoh]=V1_Amy_lfpconherency(destination,win,savelabel)
% clearvars -except destination amylfpMKII v1lfp
% savelabel=1;
% win=[301 1300];
% destination='i:\160720';
%%
tic
load([destination '\V1CONblockunitLFP.mat'])
load([destination '\Amy_unit_block_conLFP.mat'])
load([destination '\ORIvalidid.mat'])
load([destination '\XYvalidid.mat'])
load([destination '\Blksortedcell.mat'])
params.Fs = 1000; % sampling frequency
params.tapers = [2 3]; % taper parameters
params.trialave = 1;
params.fpass = [4 100];
params.pad = 0;
params.err =[2 0.05]; 

%% reform AMY lfp
celldata=cellfun(@(x) merge_block(x,win),amyunitblockconLFP,'UniformOutput',false);
unfoldeddata=cell2mat(celldata(:));
zzz=mat2cell(unfoldeddata,size(unfoldeddata,1),repmat(size(unfoldeddata,2)/numel(Blksortedcell),1,numel(Blksortedcell)));
ccc=cellfun(@(x) mat2cell(x,repmat(win(2)-win(1)+1,24,1),size(x,2)),zzz,'UniformOutput',false);
amylfpMKII=cellfun(@(x,y) x(y(:)),ccc,Blksortedcell,'UniformOutput',false);
%% V1 lfp
ccc=cellfun(@(x) merge_block(x,win),V1CONblockunitLFP,'UniformOutput',false);
Nv1validcell=sum(sum(0~=cellfun(@length,V1CONblockunitLFP{1})));
v1lfp=cellfun(@(x) mat2cell(x,size(x,1),repmat(size(x,2)/Nv1validcell,1,Nv1validcell)),ccc,'UniformOutput',false);
%% compute
LFPCoh=cellfun(@(x,y) block_coh(x,y,params),amylfpMKII',v1lfp,'UniformOutput',false);
if savelabel
    save([destination '\LFPCoh' num2str(win(1)-300) '-' num2str(win(2)-300) 'ms'],'LFPCoh','-v7.3');
end
toc
%%
function b=merge_block(a,win)
id=cellfun('length',a);
a(id==0)=[];
a=a(:);
bb=cellfun(@(x) cell2mat(x'),a,'UniformOutput',false);
b=cell2mat(bb');
b=b(win(1):win(2),:);

function C=block_coh(amy,v1,params)
C=cellfun(@(x) amy_coh(x,v1,params),amy,'UniformOutput',false);

function Coh=amy_coh(amy,v1,params)
Coh=struct;
[c1,phi1,c2,phi2,cerr1,cerr2,f1,f2]=cellfun(@(x) amy_v1_coh(amy,x,params),v1,'UniformOutput',false);
Coh.C1=cell2mat(c1);
Coh.Phi1=cell2mat(phi1);
Coh.C2=cell2mat(c2);
Coh.Phi2=cell2mat(phi2);
Coh.Cerr1=cell2mat(cerr1);
Coh.Cerr2=cell2mat(cerr2);
Coh.f1=f1;
Coh.f2=f2;

function [c1,phi1,c2,phi2,cerr1,cerr2,f1,f2]=amy_v1_coh(amy,v1,params)
neuidx=randperm(90,30)+30;
[c1,phi1,~,~,~,f1,~,~,cerr1]=coherencyc(double(amy(:,1:30)),double(v1(:,1:30)),params);
[c2,phi2,~,~,~,f2,~,~,cerr2]=coherencyc(double(amy(:,neuidx)),double(v1(:,neuidx)),params);