function solutions=randomsamplingBo_reduced(model, goodRxns, originalRxns, nSamples,supressErrors)
if nargin<4
    nSamples=1000;
end
if nargin<5
    supressErrors=false;
end

nRxns=2; %Number of reactions in the objective function in each iteration


%Reserve space for a solution matrix
sols=zeros(numel(model.rxns),nSamples);

%Main loop
counter=1;
badSolutions=0;
display(numel(goodRxns));
bms_rnx=find(strcmp('RCR99999',model.rxns));
rng(2704607)
while counter<=nSamples
   rxns=randsample(numel(goodRxns),nRxns);
   model.c=zeros(numel(model.rxns),1);
   multipliers=randsample([-1 1],nRxns,true);
   multipliers(model.rev(goodRxns(rxns))==0)=1;
   model.c([goodRxns(rxns(1));bms_rnx])=rand(nRxns,1).*multipliers;
   sol=solveLP(model);
   if any(sol.x)
       if abs(sol.f)>10^-8
            sols(:,counter)=sol.x;
            counter=counter+1;
            badSolutions=0;
       else
            badSolutions=badSolutions+1;
            %If it only finds bad solutions then throw an error.
            if badSolutions==50 && supressErrors==false
                dispEM('The program is having problems finding non-zero solutions that are not involved in loops. Review the constraints on your model. Set supressErrors to true to ignore this error'); 
            end
       end
   end
end

%Map to original model
[crap I]=ismember(model.rxns,originalRxns);
solutions=zeros(numel(originalRxns),nSamples);
solutions(I,:)=sols;
solutions=sparse(solutions);
end

%To use instead of the normal Matlab randsample function. This is in order
%to not depend on the Matlab statistical toolbox.
function I=randsample(n,k,replacement)
    if nargin<3
        replacement=false;
    end
    %n can be a integer, which leads to I being sampled from 1:n, or it can
    %be a population to sample from.
    if numel(n)==1 && isnumeric(n)
       n=1:n;
    end
    %Loop and get random numbers until the list is unique. This is only a
    %good option is the number of samples is small compared to the
    %population. There are several checks that should be made here, for
    %example regarding size and that the number of samples is <=population
    %size if replacement==false. This is not the case in randomSampling, so
    %such checks are ignored
    while true
        J=randi(numel(n),[k,1]);
        if replacement==true || numel(J)==numel(unique(J))
            I=n(J);
            break;
        end
    end
    I=I(:);
end
