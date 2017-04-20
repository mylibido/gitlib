function [S1merge,S2merge]=amyspec()
win=[401 800];
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
params.Fs = 1000; % sampling frequency
params.tapers = [2 3]; % taper parameters
params.trialave = 1;
params.fpass = [4 100];
params.pad = 0;
params.err =[2 0.05]; 
S1merge=[];
S2merge=[];
for i=1:numel(folderlist)
    tic
    load([folderlist{i} '\Amy_unit_block_conLFP.mat']);
    load([folderlist{i} '\Blksortedcell.mat']);
    celllfp=cellfun(@unfold_lfp,amyunitblockconLFP,'UniformOutput',false);
    unfoldlfp=cell2mat(celllfp(:));
    b=mat2cell(unfoldlfp,size(unfoldlfp,1),repmat(size(unfoldlfp,2)/numel(Blksortedcell),1,numel(Blksortedcell)));
    c=cellfun(@(x) mat2cell(x,repmat(size(amyunitblockconLFP{1,1}{1}{1},1),24,1),size(x,2)),b,'UniformOutput',false);
    amylfpMKII=cellfun(@(x,y) x(y(:)),c,Blksortedcell,'UniformOutput',false);
    spec=cellfun(@(x) merge_block(x,win,params),amylfpMKII,'UniformOutput',false);
   
    s1=cellfun(@(x) cell2mat(x.S1'),spec,'UniformOutput',false);
    S1=cell2mat(s1);
    s2=cellfun(@(x) cell2mat(x.S2'),spec,'UniformOutput',false);
    S2=cell2mat(s2);
    f=spec{1}.f{1};
    S1merge=[S1merge S1];
    S2merge=[S2merge S2];
    toc
end
figure
plot(f,mean(S1merge,2),'r')
hold on
plot(f,mean(S2merge,2),'g')
figure
semilogy(f,mean(S1merge,2),'r')
hold on
semilogy(f,mean(S2merge,2),'g')
sprintf('done')
function celllfp=unfold_lfp(a)
aa=cellfun(@(x) cell2mat(x'),a,'UniformOutput',false);
celllfp=cell2mat(aa');

function spec=merge_block(a,win,params)
spec=struct;
[spec.S1,spec.S2,spec.f]=cellfun(@(x) cellspec(x,win,params),a,'UniformOutput',false);
function [S1,S2,f]=cellspec(a,win,params)
neuidx=randperm(90,30)+30;
rmlfp=rmlinesc(double(a),params,[],'n',50);
[S1,f]=mtspectrumc(rmlfp(win(1):win(2),1:30),params);
[S2,~]=mtspectrumc(rmlfp(win(1):win(2),neuidx),params);


 
    