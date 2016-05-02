function [goodRxns]=combineGoodRxns(goodRxns1, goodRxns2)
B=find(goodRxns2);
A=find(goodRxns1);
d1=setdiff(A,B);
d2=setdiff(B,A);
A=[A;d2];
B=[B;d1];
A=sort(A);
B=sort(B);
if A==B
    goodRxns=A;
end
if A~=B
    dispEM('Error in good Rxns'); 
    goodRxns=[];
end
end