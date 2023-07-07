%script to generate all images by loading matlab variables and plotting and
%saving them appropriately. To be called after making changes to
%plotFigsOsher.m or plotFigsAA.
fileNames=["barbara","cameraman","mandril","disc_square"]; 

snrResults=cell(11,length(fileNames));
snrKstarResults=snrResults;
rmseResults=snrResults;
rmseKstarResults=rmseResults;

% To create and save zoomed/detail images
zoomBoxes=[[1,110,402,512];[40,110,100,170];[310,510,230,430];[30,130,120,220]];

%%
%AA images
idx=1;
for f=fileNames
    file = char(f);
    snrResults{1,idx}=file;
    snrKstarResults{1,idx}=file;
    rmseResults{1,idx}=file;
    rmseKstarResults{1,idx}=file;
    
    %load AA recoveries
    vars=dir(['./AA_blur/', file,'_noise_5/', '*.mat']); %get all .mat files
    load([vars.folder,'/',vars.name]); %load
    
    %To generate table of results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [xk_f_norm2,rmse_final,stopCrit,snr]= metricsAA(F_orig,F_data,squeeze(xkArray),T,numScales,tightFlag);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find min error and stopping criteria, original RMSE
    %Min error and index of that error. Divide by (m*n) to get RMSE
    [minVal,mink]=min(rmse_final);
    %k_star = max_k D(F_data,Txk)^2/ D(F_data,Tu)^2 \geq tau, with tau>1.
    startIdx = find(stopCrit>=1,1,'first');
    k_star = find(stopCrit(startIdx:end)<1,1,'first');
    if ~isempty(k_star)&&k_star>1
        k_star=k_star-1;
    end
    
    snrResults{2,idx}=snr(mink);
    snrKstarResults{2,idx}=snr(k_star);
    rmseResults{2,idx}=rmse_final(mink);
    rmseKstarResults{2,idx}=rmse_final(k_star);
    %end table generating code
    
    
%     %For plotting/saving figs
%     plotFigsAA(F_orig, F_data, xkArray,params,T,filePrefix,figPrefix,saveFlag,tightFlag);
    %For plotting zoomed/detailed fig crops
    zoomBox=zoomBoxes(idx,:);
    plotZoomImagesAA(F_orig, F_data, xkArray,zoomBox,T,filePrefix,figPrefix,saveFlag,tightFlag)
    
    close all; 
    clear F_orig F_data xkArray params T filePrefix figPrefix saveFlag tightFlag;

    for scheme = ["refined","tight"]
        %load AA recoveries
        if scheme == "refined"
            vars=dir(['./AA_blur_oneInit/', file,'_noise_',char(scheme),'_5/', '*.mat']); %get all .mat files
        else
            vars=dir(['./AA_blur/', file,'_noise_',char(scheme),'_5/', '*.mat']); %get all .mat files
        end
        load([vars.folder,'/',vars.name]); %load
        
        
        %for storing results in matrix for generating table or resulsts
        %%%%%%%%%%%%%%%%%%%%%%%%%
        [xk_f_norm2,rmse_final,stopCrit,snr]= metricsAA(F_orig,F_data,squeeze(xkArray),T,numScales,tightFlag);
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find min error and stopping criteria, original RMSE
        %Min error and index of that error. Divide by (m*n) to get RMSE
        [minVal,mink]=min(rmse_final);
        %k_star = max_k D(F_data,Txk)^2/ D(F_data,Tu)^2 \geq tau, with tau>1.
        startIdx = find(stopCrit>=1,1,'first');
        k_star = find(stopCrit(startIdx:end)<1,1,'first');
        if ~isempty(k_star)&&k_star>1
            k_star=k_star-1;
        end
        if scheme=="tight"
            snrResults{3,idx}=snr(mink);
            snrKstarResults{3,idx}=snr(k_star);
            rmseResults{3,idx}=rmse_final(mink);
            rmseKstarResults{3,idx}=rmse_final(k_star);
        elseif scheme=="refined"
            snrResults{4,idx}=snr(mink);
            snrKstarResults{4,idx}=snr(k_star);
            rmseResults{4,idx}=rmse_final(mink);
            rmseKstarResults{4,idx}=rmse_final(k_star);
        end
       %end table storing code
        
        
%         %uncomment to generate figs
%         plotFigsAA(F_orig, F_data, xkArray,params,T,filePrefix,figPrefix,saveFlag,tightFlag);
        
        %uncomment to generate zoomed images
        zoomBox=zoomBoxes(idx,:);
        plotZoomImagesAA(F_orig, F_data, xkArray,zoomBox,T,filePrefix,figPrefix,saveFlag,tightFlag)
    
        close all;    
        clear F_orig F_data xkArray params T filePrefix figPrefix saveFlag tightFlag;

    end
    idx=idx+1;
end
%%
%ARO images
idx=1;
for f=fileNames
    file = char(f);
    %load recoveries
    vars=dir(['./ARO_blur_optInit/', file,'_noise_additive/', '*.mat']); %get all .mat files
    load([vars.folder,'/',vars.name]); %load
    snrResults{1,idx}=file;
    snrKstarResults{1,idx}=file;
    rmseResults{1,idx}=file;
    rmseKstarResults{1,idx}=file;
    
    [xk_f_norm2,rmse_final,stopCrit,snr]= metricsARO(F_orig,F_data,squeeze(xkArray),T,numScales,tightFlag);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find min error and stopping criteria, original RMSE
    %Min error and index of that error. Divide by (m*n) to get RMSE
    [minVal,mink]=min(rmse_final);
    %k_star = max_k D(F_data,Txk)^2/ D(F_data,Tu)^2 \geq tau, with tau>1.
    startIdx = find(stopCrit>=1,1,'first');
    k_star = find(stopCrit(startIdx:end)<1,1,'first');
    if ~isempty(k_star)&&k_star>1
        k_star=k_star-1;
    end
    
    snrResults{5,idx}=snr(mink);
    snrKstarResults{5,idx}=snr(k_star);
    rmseResults{5,idx}=rmse_final(mink);
    rmseKstarResults{5,idx}=rmse_final(k_star);
    
    %For plotting and saving images
    %plotFigsARO(F_orig, F_data, xkArray,params,T,filePrefix,figPrefix,saveFlag,tightFlag) 
    
    %for plotting zoomed detail images
    zoomBox=zoomBoxes(idx,:);
    plotZoomImagesARO(F_orig, F_data, xkArray,zoomBox,T,filePrefix,figPrefix,saveFlag,tightFlag)
    
    close all; 
    clear F_orig F_data xkArray params T filePrefix figPrefix saveFlag tightFlag;
    for scheme = ["tight"]
        %load ARO recoveries
        vars=dir(['./ARO_blur_optInit/', file,'_noise_',char(scheme),'/', '*.mat']); %get all .mat files
        load([vars.folder,'/',vars.name]); %load
        
        %for storing results in matrix for generating table or resulsts
        %%%%%%%%%%%%%%%%%%%%%%%%%
        [xk_f_norm2,rmse_final,stopCrit,snr]= metricsARO(F_orig,F_data,squeeze(xkArray),T,numScales,tightFlag);
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find min error and stopping criteria, original RMSE
        %Min error and index of that error. Divide by (m*n) to get RMSE
        [minVal,mink]=min(rmse_final);
        %k_star = max_k D(F_data,Txk)^2/ D(F_data,Tu)^2 \geq tau, with tau>1.
        startIdx = find(stopCrit>=1,1,'first');
        k_star = find(stopCrit(startIdx:end)<1,1,'first');
        if ~isempty(k_star)&&k_star>1
            k_star=k_star-1;
        end
        if scheme=="tight"
            snrResults{6,idx}=snr(mink);
            snrKstarResults{6,idx}=snr(k_star);
            rmseResults{6,idx}=rmse_final(mink);
            rmseKstarResults{6,idx}=rmse_final(k_star);
        end
        %For plotting and saving images
        %plotFigsARO(F_orig, F_data, xkArray,params,T,filePrefix,figPrefix,saveFlag,tightFlag) 
    
        %for plotting zoomed detail images
        zoomBox=zoomBoxes(idx,:);
        plotZoomImagesARO(F_orig, F_data, xkArray,zoomBox,T,filePrefix,figPrefix,saveFlag,tightFlag)
        close all; 
        clear F_orig F_data xkArray params T filePrefix figPrefix saveFlag tightFlag;
    end
    idx=idx+1;
end
%%
%TNV images
idx=1;
for f=fileNames
    file = char(f);
    %load recoveries
    vars=dir(['./TNV_blur/', file,'/', '*.mat']); %get all .mat files
    load([vars.folder,'/',vars.name]); %load
    snrResults{1,idx}=file;
    snrKstarResults{1,idx}=file;
    rmseResults{1,idx}=file;
    rmseKstarResults{1,idx}=file;
    
    [xk_f_norm2,rmse_final,stopCrit,snr]= metricsTNV(F_orig,F_data,squeeze(xkArray),T,numScales,tightFlag);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find min error and stopping criteria, original RMSE
    %Min error and index of that error. Divide by (m*n) to get RMSE
    [minVal,mink]=min(rmse_final);
    %k_star = max_k D(F_data,Txk)^2/ D(F_data,Tu)^2 \geq tau, with tau>1.
    startIdx = find(stopCrit>=1,1,'first');
    k_star = find(stopCrit(startIdx:end)<1,1,'first');
    if ~isempty(k_star)&&k_star>1
        k_star=k_star-1;
    end

    snrResults{7,idx}=snr(mink);
    %snrKstarResults{10,idx}=snr(k_star);
    rmseResults{7,idx}=rmse_final(mink);
    %rmseKstarResults{10,idx}=rmse_final(k_star);
    
    %For plotting and saving images
    %plotFigsTNV(F_orig, F_data, xkArray,params,T,filePrefix,figPrefix,saveFlag,tightFlag) 
    
    %for plotting zoomed detail images. plotZoomImagesARO is okay for TNV
    %since the only difference would be in k_star computations. And TNV has
    %no meaningful kstar, so we don't consider them. 
    zoomBox=zoomBoxes(idx,:);
    plotZoomImagesARO(F_orig, F_data, xkArray,zoomBox,T,filePrefix,figPrefix,saveFlag,tightFlag)
    
    close all; 
    clear F_orig F_data xkArray params T filePrefix figPrefix saveFlag tightFlag;
    for scheme = ["tight"]
        %load ARO recoveries
        vars=dir(['./TNV_blur/', file,'_',char(scheme),'/', '*.mat']); %get all .mat files
        load([vars.folder,'/',vars.name]); %load
        
        %for storing results in matrix for generating table or resulsts
        %%%%%%%%%%%%%%%%%%%%%%%%%
        [xk_f_norm2,rmse_final,stopCrit,snr]= metricsTNV(F_orig,F_data,squeeze(xkArray),T,numScales,tightFlag);
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find min error and stopping criteria, original RMSE
        %Min error and index of that error. Divide by (m*n) to get RMSE
        [minVal,mink]=min(rmse_final);
        %k_star = max_k D(F_data,Txk)^2/ D(F_data,Tu)^2 \geq tau, with tau>1.
        startIdx = find(stopCrit>=1,1,'first');
        k_star = find(stopCrit(startIdx:end)<1,1,'first');
        if ~isempty(k_star)&&k_star>1
            k_star=k_star-1;
        end
        if scheme=="tight"
            snrResults{8,idx}=snr(mink);
            %snrKstarResults{11,idx}=snr(k_star);
            rmseResults{8,idx}=rmse_final(mink);
            %rmseKstarResults{11,idx}=rmse_final(k_star);
        end
        %For plotting and saving images
    %plotFigsTNV(F_orig, F_data, xkArray,params,T,filePrefix,figPrefix,saveFlag,tightFlag) 
    
    %for plotting zoomed detail images. plotZoomImagesARO is okay for TNV
    %since the only difference would be in k_star computations. And TNV has
    %no meaningful kstar, so we don't consider them. 
    zoomBox=zoomBoxes(idx,:);
    plotZoomImagesARO(F_orig, F_data, xkArray,zoomBox,T,filePrefix,figPrefix,saveFlag,tightFlag)
    
    close all; 
    clear F_orig F_data xkArray params T filePrefix figPrefix saveFlag tightFlag;
    end
    idx=idx+1;
end
%%
%create tables and store results as csv files
tab1=cell2table(snrResults(2:end,:),'VariableNames',snrResults(1,:));
tab2=cell2table(snrKstarResults(2:end,:),'VariableNames',snrKstarResults(1,:));
tab3=cell2table(rmseResults(2:end,:),'VariableNames',rmseResults(1,:));
tab4=cell2table(rmseKstarResults(2:end,:),'VariableNames',rmseKstarResults(1,:));
writetable(tab1,'snrResults_blur.csv')
writetable(tab2,'snrKstarResults_blur.csv')
writetable(tab3,'rmseResults_blur.csv')
writetable(tab4,'rmseKstarResults_blur.csv')


