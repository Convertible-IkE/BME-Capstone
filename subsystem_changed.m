function [fisher_result_for,fisher_result_back] = subsystem_changed(refmodel,mag,dir)
%subsystem Summary of this function goes here
% usage: [fisher_result_1,fisher_result_2] = subsystem(refmodel,minflux1,maxflux1,minflux2,maxflux2)
% by Issac Li 12/1/2014
%mag_down=find(mag<-1e-12);
mag_cutoff=prctile(mag,90);
mag_up=find(mag>mag_cutoff);
dir_cutoff1=prctile(dir,97);
dir_cutoff2=prctile(dir,3);

dir_pos=find(dir>dir_cutoff1);
dir_neg=find(dir<dir_cutoff2);

forward=intersect(mag_up,dir_pos);
backward=intersect(mag_up,dir_neg);
subsys_for=refmodel.subSystems(forward);
subsys_back=refmodel.subSystems(backward);
%%dir_cutoff1=prctile(dir,90);

RXNS_for=refmodel.rxnNames(forward);
RXNS_back=refmodel.rxnNames(backward);
pair1=[RXNS_for,subsys_for];
pair2=[RXNS_back,subsys_back];

list_1=unique(subsys_for);
list_2=unique(subsys_back);
%%
all_sub=unique(refmodel.subSystems);
C=cell(length(all_sub),1);
for i = 1:length(C);
C{i}=sum(ismember(refmodel.subSystems,all_sub(i)));
end
%%
UN=cell(length(list_2),1);
for i = 1:length(UN);
UN{i}=sum(ismember(subsys_back,list_2(i)));
end
N=cell(length(list_1),1);
for i = 1:length(N);
N{i}=sum(ismember(subsys_for,list_1(i)));
end
%% construct table on single subsystem S_N in first model
x=cell(length(list_1),1);

for i = 1:length(x);
S_N= list_1(i);
S_N_all=C{ismember(all_sub,S_N)};
S_N_req=N{i};  
x{i} = table([S_N_req;length(subsys_for)-S_N_req],[S_N_all;length(refmodel.subSystems)-S_N_all],'VariableNames',{'Forward_Up','All'},'RowNames',{S_N{1},'Not A'});
end
%% find p-values
h=zeros(length(x),1);
p=zeros(length(x),1);
stats=cell(length(x),1);
rxns_1=cell(length(x),1);
for i = 1:length(x);
[h(i) p(i) stats{i}]=fishertest(x{i,1});
selected_pairs=pair1(strcmp(x{i,1}.Properties.RowNames{1},pair1(:,2)));
rxn_ind=cellfun(@(x) find(strcmp(x,refmodel.rxnNames)) ,selected_pairs, 'UniformOutput', false);
rxn_ind=cell2mat(rxn_ind);
rxns_1{i}=horzcat(selected_pairs,num2cell(dir(rxn_ind)),num2cell(mag(rxn_ind)));
end
%%
fisher_result_for.subsystem=list_1;
fisher_result_for.H=h;
fisher_result_for.p_val=p;
%FDR=mafdr(p);
%fisher_result_back.FDR=FDR;
fisher_result_for.stats=stats;
fisher_result_for.rxns=rxns_1;
%% construct table on single subsystem S_N in second model
y=cell(length(list_2),1);
for i = 1:length(y);
S_N= list_2(i);
S_N_all=C{ismember(all_sub,S_N)};
S_N_req=UN{i};  
y{i} = table([S_N_req;length(subsys_back)-S_N_req],[S_N_all;length(refmodel.subSystems)-S_N_all],'VariableNames',{'Backward_Up','All'},'RowNames',{S_N{1},'Not A'});
end
%% find p-values
h_2=zeros(length(y),1);
p_2=zeros(length(y),1);
stats_2=cell(length(y),1);
rxns_2=cell(length(y),1);
for i = 1:length(y);
[h_2(i) p_2(i) stats_2{i}]=fishertest(y{i,1});
selected_pairs2=pair2(strcmp(y{i,1}.Properties.RowNames{1},pair2(:,2)));
rxn_ind=cellfun(@(x) find(strcmp(x,refmodel.rxnNames)) ,selected_pairs2, 'UniformOutput', false);
rxn_ind=cell2mat(rxn_ind);
rxns_2{i}=horzcat(selected_pairs2,num2cell(dir(rxn_ind)),num2cell(mag(rxn_ind)));
end
%%
fisher_result_back.subsystem=list_2;
fisher_result_back.H=h_2;
fisher_result_back.p_val=p_2;
%FDR2=mafdr(p_2);
%fisher_result_back.FDR=FDR2;
fisher_result_back.stats=stats_2;
fisher_result_back.rxns=rxns_2;


end

