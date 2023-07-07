%Removing multiplicative noise: f=u*eta, where f is degraded image, and u
% is the original image. 
% model: minimize int[lambda*( fexp(-w) + w)] + TV(w) over w. Recover
% u=exp(w). 
% Multiscale: u = u0*u1*...*uk*..., decomposition of u. w also decomposes into
% addition: w=w0+w1+...+wk+..., where wi=log(ui). We refer to the partial
% sum xk=w0+...+wk.



clear all
%for saving and retrieving images 
%for saving
folder_path="Test_Images_plus1/"; %read images with no zero values
fileNames=["barbara","cameraman","pollen","mandril","circles","geometry","disc_square"]; 
images=["barbara.png","cameraman.tif","pollen.tif","mandril_gray.tif","circles.tif","geometry.tif","disc_square.png"];
imagesPNG=["barbara.png","cameraman.png","pollen.png","mandril.png","circles.png","geometry.png","disc_square.png"];

noiseImages=["barbara_noise_02.png","cameraman_noise_02.png",...
    "pollen_noise_02.png","mandril_noise_02.png","circles_noise_02.png",...
    "geometry_noise_02.png"];
noiseImages04=["barbara_noise_04.png","cameraman_noise_04.png",...
    "pollen_noise_04.png","mandril_noise_04.png","circles_noise_04.png",...
    "geometry_noise_04.png"];%for standard deviation 0.4
%run for each image
for j=[1,2,4,7]
    close all;
    %filenames for saving
    filePrefix="SO_constInit/"+fileNames(j)+"_noise_tight/";
    figPrefix=fileNames(j)+"_";
    
    %read in image and noisy image
    F_orig=imread(char(folder_path+imagesPNG(j))); 
    F_orig=double(F_orig);
%     F_data=imread(char(folder_path+noiseImages(j)));
%     F_data=double(F_data);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %setup parameters
    [n,m]=size(F_orig);
    numScales=10;
    %algo parameters
    maxIters=1000; %time iterations in solving for wk
    dt=0.01; %0.025; %timestep
    epsilon= 0.01; %for regularizing TV
    lambda0=0.01; %intial lambda
    alp0=1;
    q=3; %for update ratio for lambda: lambda_k = lambda0*q^k;
    params=[maxIters, dt, epsilon, lambda0,q,alp0];%to pass to plotting function
    tightFlag=[1,alp0];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %Form noisy image: 
    % %%% Gamma noise %%%
    rng(10); %set random number seed so we have consistent noise across runs
    a=25; %gamma noise with mean 1, standard deviation 0.2. 
    GamNoise=gamrnd(a,1/a,size(F_orig));
    
    F_data=F_orig.*GamNoise; %multiply noise into blurred image
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Storage Arrays:
    wkArray=zeros([[m n 1], numScales]);
    xkArray=zeros([[m n 1], numScales]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Run decomposition
    xk=zeros(size(F_data));
    lambda=lambda0;
    alp=alp0;
    for k=1:numScales    
        %get decomposed piece wk. 
        wk=OsherTight(F_data,xk,dt,lambda,alp, epsilon, maxIters);
        %update xk and lambda_k
        xk=wk+xk;
        lambda=lambda* q; %alternatively, use *qk for adaptive lambda
        alp = 1/(k)^(3/2);
        %Store images:
        wkArray(:,:,1,k)=wk; %image piece, single scale
        xkArray(:,:,1,k)=exp(xk); %updated multiscale image

    %     % For adaptive lambda updates, specifies how to choose qk
    %     val=D_f_data_Txk_D_f_Data_f_orig(k); %take last bregman ratio
    %     val0=D_f_data_Txk_D_f_Data_f_orig(1);%initial bregman ratio
    %     qk = 1.5/(1+5*exp(-(val-1)))^10+1;  %sigmoid qk update
    %     %qk=1+log(val); %log qk update
    %     %qk=2.5/(val0-1)*(val-1)+1; %linear fit between (1,1), (val0,2)
    %     qk=min(qk,2.0); %choose bounds of qk
    %     qk=max(qk,1.05);
    end

    %Plot and save image
    saveFlag=1;
    if saveFlag==1
        mkdir(char(filePrefix));
        save(filePrefix+figPrefix+"vars",'F_orig', 'F_data', 'xkArray','params','filePrefix','figPrefix','saveFlag','tightFlag', 'numScales')
    end
    plotFigsOsher(F_orig, F_data, xkArray,params,filePrefix,figPrefix,saveFlag, tightFlag);
    %[xk_f_norm2,rmse_final,stopCrit,snr]= metrics(F_orig+1,F_data+1,squeeze(xkArray)+1,numScales,tightFlag);

end