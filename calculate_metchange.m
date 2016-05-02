%% Read Expression Files
 direct=dir('C:\myFiles\Capstone\Data\tox_hep_rno');
 names=struct2cell(direct)';
 names=names(3:end,1);
 %% Locate Carbon Tetrachloride and Calculate  MetChange scores
 load( 'C:\myFiles\Capstone\Model\cobra_models.mat');
 addpath 'C:\myFiles\Capstone\Data\tox_hep_rno'
 addpath 'C:\myFiles\Capstone\matlab'
 ctci=~cellfun(@isempty,regexp(names,'.*carbontetrachloride.*'));
ctc=names(ctci);
for i=1:length(ctc)
     display(['running: ',num2str(i)]);    
     expressionTable=readtable(ctc{i},'Delimiter','\t');
     if iscell(expressionTable.pval(1))
          A=(cellfun(@str2num,expressionTable.pval,'UniformOutput',false));
         expressionTable((cellfun(@isempty,A)),:)=[];
         expressionTable.pval=cell2mat(cellfun(@str2num,expressionTable.pval,'UniformOutput',false));
     end
       if iscell(expressionTable.logfc(1))
         expressionTable.logfc=cell2mat(cellfun(@str2num,expressionTable.logfc,'UniformOutput',false));
     end
     filename=regexprep(ctc{i},'.txt','');
    [objDM] = MetCHANGE(rno_cobra,expressionTable.ftest_pval,'gurobi5');
    eval([filename,'=objDM;']);
end

cd 'C:\myFiles\Capstone\Matlab\'
save metchange_ccl4.mat rno_hep_t1_one_carbontetrachloride_d1_gene_efit_robust rno_hep_t1_one_carbontetrachloride_d2_gene_efit_robust rno_hep_t2_one_carbontetrachloride_d1_gene_efit_robust ...
    rno_hep_t2_one_carbontetrachloride_d2_gene_efit_robust rno_hep_t2_one_carbontetrachloride_d3_gene_efit_robust ...
    rno_hep_t3_one_carbontetrachloride_d1_gene_efit_robust rno_hep_t3_one_carbontetrachloride_d2_gene_efit_robust