clear
clc
warning off;
for index=1
DataName{1} = 'YALE';

lambdaset = 2.^[-15:1:15];
betaset = 10.^[-15:1:15];
% gammaset = 10.^[-5:2:5];

DataHyperParam{10} = [
    6,5,5;
    ];

count = 1;
for i=1:length(lambdaset)
    for j=1:length(betaset)
        allsort(count,:) = [i,j];
        count = count + 1;
    end
end

DataHyperParam{index} = allsort;

;

path = './';
addpath(genpath(path));
k=0;

    k=k+1;
    dataName = DataName{index} %%% flower17; flower102; proteinFold,caltech101_mit,UCI_DIGIT,ccv
    dataHyperParam = DataHyperParam{index};
    load([path,'datasets/',dataName,'_Kmatrix'],'KH','Y');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numclass = length(unique(Y));
    numker = size(KH,3);
    num = size(KH,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    KH = kcenter(KH);
    KH = knorm(KH);
    K0 = zeros(num,num); 
    qnorm = 2;
    opt.disp = 0;
    tic
    
    for param=1:size(dataHyperParam,1)
        it = dataHyperParam(param,1);
        ij = dataHyperParam(param,2);
        disp(['p1=',num2str(lambdaset(it)), '  p2=',num2str(betaset(ij))]);
        tic;
        
        [H_normalized9, iter, obj] = PNOF(KH,numclass,lambdaset(it),betaset(ij),Y,numclass);
        res9 = ClusteringMeasure_new(Y, H_normalized9);
        time(k)=toc;
        accval9(it,ij) = res9(1,1);      
        nmival9(it,ij)= res9(1,2);      
        purval9(it,ij) = res9(1,3);
    end
    res = [max(max(max(accval9))); max(max(max(nmival9)));max(max(max(purval9)))];
end


