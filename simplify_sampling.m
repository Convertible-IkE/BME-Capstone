function [model_up, goodRxns, originalRxns]=simplify_sampling(model)
%First check that the model is feasible given the constraints
sol=solveLP(model);
if isempty(sol.x)
   dispEM('The model has no feasible solution'); 
end 

%Simplify the model to speed stuff up a little. Keep original mapping
originalRxns=model.rxns;

%Then change the bounds to +/- Inf. This is needed in order to not have
%loops in the solutions
if nargin<2
   model.ub(model.ub==max(model.ub))=Inf;
   model.lb(model.lb==min(model.lb))=-Inf;
end

%Reactions which can be involved in loops should not be optimized
%for. Check which reactions reach an arbitary high upper bound
goodRxns=true(numel(model.rxns),1);
for i=1:numel(model.rxns)
    if goodRxns(i)==true
        testModel=setParam(model,'eq',model.rxns(i),1000);
        sol=solveLP(testModel);
        if ~isempty(sol.f)
            goodRxns(abs(sol.x)>999)=false;
        else
            %If the reaction is reversible, also check in that direction
            if model.rev(i)
                testModel=setParam(model,'eq',model.rxns(i),-1000);
                sol=solveLP(testModel);
                if ~isempty(sol.f)
                    goodRxns(abs(sol.x)>999)=false;
                end
            end
        end
    end
end
model_up=model;
end