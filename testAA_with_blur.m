% Removing multiplicative blurring and noise: f_data=Tu*eta, where f_data is degraded image, and u
% is the original image. 
% model: minimize int[lambda*( log(Tuxk) + f/Tuxk] + ak lambda TV(log(uxk))+ TV(log(u)) over u. 
% Multiscale: u = u0*u1*...*uk*..., decomposition of u. w We refer to
% partial prodcut xk=u0*u1*...*uk

clear all

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


for j=[2] %loop over all images
    close all;
    %filenames for saving
    filePrefix="AA_blur/"+fileNames(j)+"_noise_5/"; %for blur+noise
    %filePrefix="AA_blur/"+fileNames(j)+"/"; %for blur only
    mkdir(char(filePrefix));
    figPrefix=fileNames(j)+"_";

    %read in image and noisy image
    F_orig=imread(char(folder_path+imagesPNG(j))); 
    F_orig=double(F_orig);
%     F_data=imread(char(folder_path+noiseImages(j)));
%     F_data=double(F_data);
    
    %setup parameters
    [n,m]=size(F_orig);
    numScales=2;
    %algo parameters
    maxIters=4000; %time iterations in solving for wk
    dt=0.01; %0.025; %timestep
    epsilon= 0.01; %for regularizing TV
    lambda0=0.01; %intial lambda
    alp0=1;
    q=2; %for update ratio for lambda: lambda_k = lambda0*q^k;
    params=[maxIters, dt, epsilon, lambda0,q,alp0];%to pass to plotting function
    tightFlag=[0,alp0];% to pass to plotting and metrics functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Form noisy image: switch T to identity for only noise. Switch GamNoise to
    %ones(size(F_orig)) for only blurring. 
    %%% Blur and Gamma noise %%%
    rng(10); %fix seed across all runs
    a=25; %gamma noise with mean 1, standard deviation 0.2. 
    GamNoise=gamrnd(a,1/a,size(F_orig)); %noise
    %GamNoise=ones(size(F_orig)); %For only deblurring, no noise

    T=fspecial('gaussian',[5 5],sqrt(2)); %blurring component/operator
    %T=fspecial('average',[1 1]); %identity, for no blur
    F_blur=imfilter(F_orig,T,'symmetric','same'); %create blurred image
    F_data=F_blur.*GamNoise; %multiply noise into blurred image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Storage Arrays:
    ukArray=zeros([[m n 1], numScales]);
    xkArray=zeros([[m n 1], numScales]);
    energies=zeros(maxIters,numScales);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Run decomposition
    xk=ones(size(F_data));
    lambda=lambda0;
    for k=1:numScales    
        %get decomposed piece uk. 
        [uk,energy]=AAlog_blur(F_data,xk,dt,lambda,T,epsilon,maxIters);
        %update xk and lambda_k
        xk=uk.*xk;
        lambda=lambda* q; %alternatively, use *qk for adaptive lambda

        %Store images:
        ukArray(:,:,1,k)=uk; %image piece, single scale
        xkArray(:,:,1,k)=xk; %updated multiscale image
        energies(1:length(energy),k) = energy;

    %     % For adaptive lambda updates, specifies how to choose qk
    %     val=D_f_data_Txk_D_f_Data_f_orig(k); %take last bregman ratio
    %     val0=D_f_data_Txk_D_f_Data_f_orig(1);%initial bregman ratio
    %     qk = 1.5/(1+5*exp(-(val-1)))^10+1;  %sigmoid qk update
    %     %qk=1+log(val); %log qk update
    %     %qk=2.5/(val0-1)*(val-1)+1; %linear fit between (1,1), (val0,2)
    %     qk=min(qk,2.0); %choose bounds of qk
    %     qk=max(qk,1.05);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot and save
    saveFlag=0; %set to 1 to save, leave 0 to just display
    if saveFlag==1
        save(filePrefix+figPrefix+"vars",'F_orig', 'F_data', 'xkArray','params','T','filePrefix','figPrefix','saveFlag','tightFlag', 'numScales')
    end
    plotFigsAA(F_orig, F_data, xkArray,params,T,filePrefix,figPrefix,saveFlag,tightFlag);
    [xk_f_norm2,rmse_final,stopCrit,snr]= metricsAA(F_orig,F_data,squeeze(xkArray),T,numScales,tightFlag);
    figure()
    for k=1:numScales
        hold on; 
        vec=energies(energies(:,k)>0,k);
        semilogy(1:length(vec),vec);
    end
    legend({'1','2','3','4','5','6','7','8','9','10','11','12','13','14'})

end