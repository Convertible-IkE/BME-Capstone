function [dir,dir_std,delta_mag,delta_std,mag]=TwoDSystemPlot(SparsematrixDrug,SparsematrixCtrl,PlotTitle)
 fullMatrixA=full(SparsematrixDrug );
 fullMatrixB=full(SparsematrixCtrl);
ctrl_mean=mean(fullMatrixB(:,:),2);
drug_mean=mean(fullMatrixA(:,:),2);
dir=(drug_mean-ctrl_mean);
ctrl_sqsum=sum(fullMatrixB.^2,2);
drug_sqsum=sum(fullMatrixA.^2,2);
delta_mag=sqrt(abs(drug_sqsum-ctrl_sqsum))/size(fullMatrixB,2);
mag=(sum(abs(fullMatrixA)+abs(fullMatrixB),2))/size(fullMatrixB,2);
dir_std=dir./mag;
delta_std=delta_mag./mag;
scatter(dir,delta_mag);
title(PlotTitle)
xlabel('Direction')
ylabel('Delta Magnitude')
figure
scatter(delta_std,dir_std);
ylabel('Normalized Direction')
xlabel('Normalized Magnitude')
title(PlotTitle)
ylim([-1,1])
end