%Removing multiplicative noise: f=Tu*eta, where f is degraded image, and u
% is the original image. 
% model: minimize int[lambda*(fexp(-ux_k)+ux_k] + TV(u) over u. 
% Multiscale: u = u0*u1*...*uk*..., decomposition of u. We refer to the partial
% produce xk=u0*...*uk.
clear all
%for saving
folder_path="Test_Images_plus1/"; %read images with no zero values
fileNames=["barbara","cameraman","pollen","mandril","circles","geometry","disc_square"]; 
images=["barbara.png","cameraman.tif","pollen.tif","mandril_gray.tif","circles.tif","geometry.tif","disc_square.png"];
imagesPNG=["barbara.png","cameraman.png","pollen.png","mandril.png","circles.png","geometry.png","disc_square.png"];

ims=[1,2,4,7];
for j=ims %1:length(images) %loop over all images
    close all;
    %filenames for saving
    filePrefix="AA_MHDM/"+fileNames(j)+"_tight/";
    figPrefix=fileNames(j)+"_";

    %read in image and noisy image
    F_orig=imread(char(folder_path+imagesPNG(j))); 
    F_orig=double(F_orig);
%     F_data=imread(char(folder_path+noiseImages(j)));
%     F_data=double(F_data);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %setup parameters
    [n,m]=size(F_orig);
    numScales=80;
    %algo parameters
    maxIters=300; %time iterations in solving for wk
    dt=0.01; %0.025; %timestep
    epsilon= 0.01; %for regularizing TV
    lambda0=0.01; %intial lambda
    alp0=1;
    q=3; %for update ratio for lambda: lambda_k = lambda0*q^k;
    params=[maxIters, dt, epsilon, lambda0,q,alp0];%to pass to plotting function
    tightFlag=[1,alp0]; %to pass to plotting and metrics functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Form noisy image: 
    %%% Gamma noise %%%
    rng(10);
    a=25; %gamma noise with mean 1, standard deviation 0.2. 
    GamNoise=gamrnd(a,1/a,size(F_orig));
    
    %T=fspecial('gaussian',[3 3],sqrt(2)); %blurring component/operator
    T=fspecial('average',[1 1]); %identity, for no blur
    F_blur=imfilter(F_orig,T,'symmetric','same'); %create blurred image
    F_data=F_blur.*GamNoise; %multiply noise into blurred image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Storage Arrays:
    ukArray=zeros([[m n 1], numScales]);
    xkArray=zeros([[m n 1], numScales]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Run decomposition
    xk=ones(size(F_data));
    lambda=lambda0;
    alp=alp0;
    for k=1:numScales    
        %get decomposed piece wk. 
        uk=AA_tight_MHDM(F_data,xk,dt,lambda,alp,T, epsilon, maxIters);
        %update xk and lambda_k
        xk=uk.*xk;
        lambda=lambda* q; %alternatively, use *qk for adaptive lambda
        alp=1/(k^(3/2));
        %Store images:
        ukArray(:,:,1,k)=uk; %image piece, single scale
        xkArray(:,:,1,k)=xk; %updated multiscale image

    %     % For adaptive lambda updates, specifies how to choose qk
    %     val=D_f_data_Txk_D_f_Data_f_orig(k); %take last bregman ratio
    %     val0=D_f_data_Txk_D_f_Data_f_orig(1);%initial bregman ratio
    %     qk = 1.5/(1+5*exp(-(val-1)))^10+1;  %sigmoid qk update
    %     %qk=1+log(val); %log qk update
    %     %qk=2.5/(val0-1)*(val-1)+1; %linear fit between (1,1), (val0,2)
    %     qk=min(qk,2.0); %choose bounds of qk
    %     qk=max(qk,1.05);
    end

    %Plot 
    saveFlag=1;
    if saveFlag==1
        mkdir(char(filePrefix));
        save(filePrefix+figPrefix+"vars",'F_orig', 'F_data', 'xkArray','params','T','filePrefix','figPrefix','saveFlag','tightFlag', 'numScales')
    end
    plotFigsAA(F_orig, F_data, xkArray,params,T,filePrefix,figPrefix,saveFlag,tightFlag)

    %to get metrics for inspection
    %[xk_f_norm2,rmse_final,stopCrit,snr]= metrics(F_orig+1,F_data+1,squeeze(xkArray)+1,numScales,tightFlag);
end