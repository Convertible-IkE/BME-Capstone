A=rno_hep_t3_one_carbontetrachloride_d2_gene_efit_robust;

%%
A_m1=mean(A.ctrl(:,:),2);
A_m2=mean(A.drug(:,:),2);

for i=1:length(A_m1)
   % if ((abs(A_m2(i))>5)+ (abs(A_m2(i))<10)==2)
    if (abs(A_m2(i))>500)
   figure
    subplot(2,1,1)
    hist(A.ctrl(i,:),20)
    title('Ctrl')
    subplot(2,1,2)
    hist(A.drug(i,:),20)
    title('Drug Treament')

    end
end

%%
met_rxncount=zeros(1,1);
for i=1:length(rno_cobra.mets)
met_rxncount(i,1)=sum(rno_cobra.S(i,:)~=0);
end
hist(met_rxncount)
ph = get(gca,'children');
% Determine number of histogram patches
N_patches = length(ph);
for i = 1:N_patches
      % Get patch vertices
      vn = get(ph(i),'Vertices');
      % Adjust y location
      vn(:,2) = vn(:,2) + 1;
      % Reset data
      set(ph(i),'Vertices',vn)
end
% Change scale
set(gca,'yscale','log')
title('Ctrl and Drug Treatment')

%%
for i=1:length(rno_cobra.mets)
    met_rxns(i,1)=sum(rno_cobra.S(i,:)~=0);
    met_rxns(i,2)=sum(rno_cobra.S(i,:)>0);
    met_rxns(i,3)=sum(rno_cobra.S(i,:)<0);
    met_pos_rxns{i,1}=rno_cobra.rxns(rno_cobra.S(i,:)>0);
    met_neg_rxns{i,1}=(rno_cobra.rxns(rno_cobra.S(i,:)<0));    
end
    met_rxns=full(met_rxns);
rno_cobra.rxns(rno_cobra.S(1,:)~=0)

%% Caculate Net Metabolite Fluxes and Compare between treatments
% use mean as a toy
A_m1=A_m1;
w1=diag(full(A_m1));
w2=diag(full(A_m2));
S1=full(rno_cobra.S);
S2=S1;
S1=S1*w1;
S2=S2*w2;
c_d=diag(~c);
S1=S1*c_d;
S2=S2*c_d;
mf1=sum(S1,2);
mf2=sum(S2,2);

%% flux matrix
 direct=dir('C:\myFiles\Capstone\Data\tox_hep_rno');
 names=struct2cell(direct)';
 names=names(3:end,1);
ctci=~cellfun(@isempty,regexp(names,'.*carbontetrachloride.*'));
ctc=names(ctci);
ctc=cellfun(@(x) strrep(x,'.txt','') ,ctc, 'UniformOutput',false);
cormatrix=[];
% Matrix with direct means of flux
for i=1:length(ctc)
 cormatrix=[cormatrix, eval(['mean(',ctc{i},'.drug,2)'])];
 
  cormatrix=[cormatrix, eval(['mean(',ctc{i},'.ctrl,2)'])];
end

imagesc(A)
colormap(jet);
colorbar;

% Matrix with S transform
S1=full(rno_cobra.S);
S2=S1;
mfmatrix=[];
 [ex,ut]=findExcRxns(rno_cobra);
 sec=ex-ut;
 sec_d=diag(~sec);
for i=1:length(ctc)
 A_m1=eval(['mean(',ctc{i},'.ctrl,2)']);
 A_m2=eval(['mean(',ctc{i},'.drug,2)']);
 w1=diag(full(A_m1));
w2=diag(full(A_m2));
S1=full(rno_cobra.S);
S2=S1;
S1=S1*w1;
S2=S2*w2;
S1=S1*sec_d;
S2=S2*sec_d;
mf1=sum(S1,2);
mf2=sum(S2,2);
mfmatrix=[mfmatrix,mf1,mf2];
end

%%
S1=full(rno_cobra.S);
S2=S1;
mfmatrix=[];
 [ex,ut]=findExcRxns(rno_cobra);
 ctc_t3_d2_drug=full(rno_hep_t3_one_carbontetrachloride_d2_gene_efit_robust.drug);
 ctc_t2_d3_drug=full(rno_hep_t2_one_carbontetrachloride_d3_gene_efit_robust.drug);
 ctc_t3_d2_ctrl=full(rno_hep_t2_one_carbontetrachloride_d3_gene_efit_robust.ctrl);
 ctc_t2_d3_ctrl=full(rno_hep_t2_one_carbontetrachloride_d3_gene_efit_robust.ctrl);
 
 ex_ctc_t3_d2_drug=ctc_t2_d3_drug(ex,:);
 ex_ctc_t3_d2_ctrl=ctc_t2_d3_ctrl(ex,:);
 %%
ctrl_mean=mean(ex_ctc_t3_d2_ctrl(:,:),2);
drug_mean=mean(ex_ctc_t3_d2_drug(:,:),2);
dir=(drug_mean-ctrl_mean);    
ctrl_sqsum=sum(ex_ctc_t3_d2_ctrl.^2,2);
drug_sqsum=sum(ex_ctc_t3_d2_drug.^2,2);
mag=sqrt(abs(drug_sqsum-ctrl_sqsum)/size(ex_ctc_t3_d2_ctrl,2));
display(rno_cobra.rxnNames(ex))

%%
ctrl_mean=mean(ctc_t3_d2_ctrl(:,:),2);
drug_mean=mean(ctc_t3_d2_drug(:,:),2);
dir=(drug_mean-ctrl_mean);    
ctrl_sqsum=sum(ctc_t3_d2_ctrl.^2,2);
drug_sqsum=sum(ctc_t3_d2_drug.^2,2);
mag=sqrt(abs(drug_sqsum-ctrl_sqsum)/size(ctc_t3_d2_ctrl,2));
display(rno_cobra.rxnNames(ex))


%%
index=(mag<6000)+(mag>1000);
rno_cobra.rxnNames(index==2)
%find(index)


subplot(2,1,1)
hist(ctc_t3_d2_ctrl(end,:))
title('Ctrl')
subplot(2,1,2)
hist(ctc_t3_d2_drug(end,:))
title('Drug')


help made