%script to generate metric plots by loading matlab variables and plotting and
%saving them appropriately. To be called after making changes to
%plotFigsOsher.m or plotFigsAA.
clear all;
fileNames=["cameraman"]; 

kstars=zeros(1,3);
kmins=zeros(1,3);
% rmses;
% snrs;
% stopCrits;
%Osher shi images
idx=1;
for f=fileNames
    file = char(f);
    %load osher shi recoveries
    vars=dir(['./additive/', file,'_noise_additive/', '*.mat']); %get all .mat files
    load([vars.folder,'/',vars.name]); %load    
    [xk_f_norm2,rmse_final,stopCrit,snr]= metrics(F_orig,F_data,squeeze(xkArray),numScales,tightFlag);
    rmses{1}=rmse_final;
    snrs{1}=snr;
    stopCrits{1}=stopCrit;
    
    for scheme = ["tight","refined"]
        %load osher shi recoveries
        vars=dir(['./additive/', file,'_noise_',char(scheme),'/', '*.mat']); %get all .mat files
        load([vars.folder,'/',vars.name]); %load
        [xk_f_norm2,rmse_final,stopCrit,snr]= metrics(F_orig,F_data,squeeze(xkArray),numScales,tightFlag);
        %
        if scheme=="tight"
            rmses{2}=rmse_final;
            snrs{2}=snr;
            stopCrits{2}=stopCrit;
        elseif scheme=="refined"
            rmses{3}=rmse_final;
            snrs{3}=snr;
            stopCrits{3}=stopCrit;
        end
    end
    
    for k=1:3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find min error and stopping criteria, original RMSE
        %Min error and index of that error. Divide by (m*n) to get RMSE
        [minVal,mink]=min(rmses{k});
        %k_star = max_k D(F_data,Txk)^2/ D(F_data,Tu)^2 \geq tau, with tau>1.
        k_star=min(find((stopCrits{k}<=1)==1));
        if ~isempty(k_star)&&k_star>1
            k_star=k_star-1;
        end
        kstars(1,k)=k_star;
        kmins(1,k)=mink;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    %Plot the RMSE and SNR, as well as the stopping criterion plot
    fig=figure('position',[100,100,1000,500]);
    p=10; % markersize
    l=1; %linesize
    set(fig,'defaultAxesColorOrder',[[0,0,0];[0,0,0]]);
    hold on;
    
    ax=gca; ax.FontSize=16;
    yyaxis left
    p1=plot(0:length(rmses{1})-1,rmses{1},'-d','LineWidth',l,'MarkerSize',p);
    p2=plot(0:length(rmses{2})-1,rmses{2},'-v','LineWidth',l,'MarkerSize',p);
    p3=plot(0:length(rmses{3})-1,rmses{3},'-+','LineWidth',l,'MarkerSize',p);
    
    yyaxis left
    plot(kstars(1)-1,rmses{1}(kstars(1)),'r-*','LineWidth',l+1,'MarkerSize',p-2)
    plot(kstars(2)-1,rmses{2}(kstars(2)),'r-*','LineWidth',l+1,'MarkerSize',p-2)
    plot(kstars(3)-1,rmses{3}(kstars(3)),'r-*','LineWidth',l+1,'MarkerSize',p-2)
    
    plot(kmins(1)-1,rmses{1}(kmins(1)),'b-s','LineWidth',l+1,'MarkerSize',p)
    plot(kmins(2)-1,rmses{2}(kmins(2)),'b-s','LineWidth',l+1,'MarkerSize',p)
    plot(kmins(3)-1,rmses{3}(kmins(3)),'b-s','LineWidth',l+1,'MarkerSize',p)

    
    xlabel('Multiscales: k','FontSize',22)
    ylabel('RMSE','FontSize',22)
    title('RMSE and SNR vs Multiscales','FontSize',22)
    
    yyaxis right
    p4=plot(0:length(snrs{1})-1,snrs{1},'r-.d','LineWidth',l,'MarkerSize',p);
    p5=plot(0:length(snrs{2})-1,snrs{2},'r-.v','LineWidth',l,'MarkerSize',p);
    p6=plot(0:length(snrs{3})-1,snrs{3},'r-.+','LineWidth',l,'MarkerSize',p)  ;
%     plot(kstars(1),snrs{1}(kstars(1)),'g-o','LineWidth',3)
%     plot(kstars(2),snrs{2}(kstars(2)),'g-o','LineWidth',3)
%     plot(kstars(3),snrs{3}(kstars(3)),'g-o','LineWidth',3)
%     
%     plot(kmins(1),snrs{1}(kmins(1)),'m-o','LineWidth',3)
%     plot(kmins(2),snrs{2}(kmins(2)),'m-o','LineWidth',3)
%     plot(kmins(3),snrs{3}(kmins(3)),'m-o','LineWidth',3)
%     
    
    ylabel('SNR','FontSize',22)
    legend([p1,p2,p3,p4,p5,p6],{'RMSE MHDM','RMSE Tight','RMSE Refined','SNR MHDM','SNR Tight','SNR Refined'},'FontSize',18)
    
 if saveFlag==1
        figName="./General_Figs/combined_metrics.png";
        saveas(gcf,figName)
        im=imread(char(figName));
        im_crop=crop_borders(im,[255],[0.01]);
        imwrite(uint8(im_crop),char(figName));
 end
  
figure('position',[100,100,700,500])
hold on
ax=gca; ax.FontSize=16;
semilogy(0:length(stopCrits{1})-1, stopCrits{1},'k-d','LineWidth',l,'MarkerSize',p)
semilogy(0:length(stopCrits{2})-1, stopCrits{2},'k-v','LineWidth',l,'MarkerSize',p)
semilogy(0:length(stopCrits{3})-1, stopCrits{3},'k-+','LineWidth',l,'MarkerSize',p)
semilogy(0:length(stopCrits{1})-1,ones(length(stopCrits{1}),1),'k--','LineWidth',2)
plot(kstars(1)-1,stopCrits{1}(kstars(1)),'r-*','LineWidth',l+1,'MarkerSize',p-2)
plot(kstars(2)-1,stopCrits{2}(kstars(2)),'r-*','LineWidth',l+1,'MarkerSize',p-2)
plot(kstars(3)-1,stopCrits{3}(kstars(3)),'r-*','LineWidth',l+1,'MarkerSize',p-2)
xlabel('Multiscales: k','FontSize',22)
title('Stopping Criteria','FontSize',22)
legend({'MHDM','Tight MHDM','Refined MHDM'},'FontSize',18)
%     if tightFlag(1)==1
%         title(['(D(F_{data},x_k)+\alpha_kTV(log(x_k)))/D(F_{data},u), k^*=',num2str(k_star)],'FontSize',22)    
%     else
%         title(['D(F_{data},x_k)/D(F_{data},u), k^*=',num2str(k_star)],'FontSize',22)
%     end

    if saveFlag==1
        figName="./General_Figs/combined_stopping.png";
        saveas(gcf,figName)
        im=imread(char(figName));
        im_crop=crop_borders(im,[255],[0.01]);
        imwrite(uint8(im_crop),char(figName));
    end
%     
    idx=idx+1;
end