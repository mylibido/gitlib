function [LFPCoh]=V1_Amy_lfpcoherency(destination,win,shufflelabel,savelabel)
% [LFPCoh]=V1_Amy_lfpconherency(destination,win,savelabel)
% shufflelabel =0;
% savelabel=[1 1 1];%  save Coh Amy V1
% win=[301 1300];
% destination='i:\160720';
%%
if numel(savelabel)~=3
    error('savelabel should be like [1 0 0]');
end 
tic
shuffleL=shufflelabel;
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
if ~exist([destination '\AmylfpMKII' num2str(win(1)-300) '-' num2str(win(2)-300) 'ms.mat'],'file')
    celldata=cellfun(@(x) merge_block(x,win,params),amyunitblockconLFP,'UniformOutput',false);
    unfoldeddata=cell2mat(celldata(:));
    zzz=mat2cell(unfoldeddata,size(unfoldeddata,1),repmat(size(unfoldeddata,2)/numel(Blksortedcell),1,numel(Blksortedcell)));
    ccc=cellfun(@(x) mat2cell(x,repmat(win(2)-win(1)+1,24,1),size(x,2)),zzz,'UniformOutput',false);
    AmylfpMKII=cellfun(@(x,y) x(y(:)),ccc,Blksortedcell,'UniformOutput',false);
else
    load([destination '\AmylfpMKII' num2str(win(1)-300) '-' num2str(win(2)-300) 'ms.mat'])
end
%% V1 lfp
if ~exist([destination '\V1lfp' num2str(win(1)-300) '-' num2str(win(2)-300) 'ms.mat'],'file')
    ccc=cellfun(@(x) merge_block(x,win,params),V1CONblockunitLFP,'UniformOutput',false);
    Nv1validcell=sum(sum(0~=cellfun(@length,V1CONblockunitLFP{1})));
    V1lfp=cellfun(@(x) mat2cell(x,size(x,1),repmat(size(x,2)/Nv1validcell,1,Nv1validcell)),ccc,'UniformOutput',false);
else
    load([destination '\V1lfp' num2str(win(1)-300) '-' num2str(win(2)-300) 'ms.mat'])
end
%% compute
LFPCoh=cellfun(@(x,y) block_coh(x,y,params,shuffleL),AmylfpMKII',V1lfp,'UniformOutput',false);
if savelabel(1)
    if shuffleL
    save([destination '\LFPCoh_shuffle' num2str(win(1)-300) '-' num2str(win(2)-300) 'ms'],'LFPCoh','-v7.3');    
    else
    save([destination '\LFPCoh' num2str(win(1)-300) '-' num2str(win(2)-300) 'ms'],'LFPCoh','-v7.3');
    end
end
if savelabel(2)
    save([destination '\AmylfpMKII' num2str(win(1)-300) '-' num2str(win(2)-300) 'ms'],'AmylfpMKII','-v7.3');
end
if savelabel(3)
    save([destination '\V1lfp' num2str(win(1)-300) '-' num2str(win(2)-300) 'ms'],'V1lfp','-v7.3');
end
toc
%%
function b=merge_block(a,win,params)
id=cellfun('length',a);
a(id==0)=[];
a=a(:);
bb=cellfun(@(x) cell2mat(x'),a,'UniformOutput',false);
bblsc=cellfun(@(x) rmlinesc(double(x),params,[],'n',50),bb,'UniformOutput',false);
b=cell2mat(bblsc');
b=b(win(1):win(2),:);

function C=block_coh(amy,v1,params,shuffleL)
C=cellfun(@(x) amy_coh(x,v1,params,shuffleL),amy,'UniformOutput',false);

function Coh=amy_coh(amy,v1,params,shuffleL)
Coh=struct;
[c1,phi1,c2,phi2,cerr1,cerr2,f1,f2]=cellfun(@(x) amy_v1_coh(amy,x,params,shuffleL),v1,'UniformOutput',false);
Coh.C1=cell2mat(c1);
Coh.Phi1=cell2mat(phi1);
Coh.C2=cell2mat(c2);
Coh.Phi2=cell2mat(phi2);
Coh.Cerr1=cell2mat(cerr1);
Coh.Cerr2=cell2mat(cerr2);
Coh.f1=f1;
Coh.f2=f2;

function [c1,phi1,c2,phi2,cerr1,cerr2,f1,f2]=amy_v1_coh(amy,v1,params,shuffleL)
neuidx=randperm(90,30)+30;
if shuffleL
    Nshuffle=20;
    randtim=1:Nshuffle;
    randidx=arrayfun(@(x) randperm(30),randtim,'UniformOutput',false);
    [c1c,phi1c,~,~,~,f1c,~,~,cerr1c]=cellfun(@(x) coherencyc(amy(:,1:30),v1(:,x),params),randidx,'UniformOutput',false);
    [c2c,phi2c,~,~,~,f2c,~,~,cerr2c]=cellfun(@(x) coherencyc(amy(:,neuidx),v1(:,neuidx(x)),params),randidx,'UniformOutput',false);
    c1=mean(cell2mat(reshape(c1c,1,1,[])),2);
    phi1=mean(cell2mat(reshape(phi1c,1,1,[])),2);
    cerr1=mean(cell2mat(reshape(cerr1c,1,1,[])),2);
    f1=f1c{1};
    c2=mean(cell2mat(reshape(c2c,1,1,[])),2);
    phi2=mean(cell2mat(reshape(phi2c,1,1,[])),2);
    cerr2=mean(cell2mat(reshape(cerr2c,1,1,[])),2);
    f2=f2c{1};
else
    [c1,phi1,~,~,~,f1,~,~,cerr1]=coherencyc(amy(:,1:30),v1(:,1:30),params);
    [c2,phi2,~,~,~,f2,~,~,cerr2]=coherencyc(amy(:,neuidx),v1(:,neuidx),params);
end