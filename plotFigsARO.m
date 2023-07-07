function plotFigsARO(F_orig, F_data, xkArray,params,T,filePrefix,figPrefix,saveFlag,tightFlag)

[m,n]=size(F_orig);
numScales = length(xkArray(1,1,1,:));
%tightFlag=[0,0]; %tightFlag(1)=indicates tight or not, tightFlag(2)= value of alp0
if length(params)==6
    paramsCell=num2cell(params);
    [maxIters, dt, epsilon, lambda0, q,alp0]=paramsCell{:};
else 
    paramsCell=num2cell(params);
    [maxIters, dt, epsilon, lambda0,q]=paramsCell{:};
end
[xk_f_norm2,rmse_final,stopCrit,snr]= metricsARO(F_orig,F_data,squeeze(xkArray),T,numScales,tightFlag);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find min error and stopping criteria, original RMSE
%Min error and index of that error. Divide by (m*n) to get RMSE
[minVal,mink]=min(rmse_final);
%k_star = max_k D(F_data,Txk)^2/ D(F_data,Tu)^2 \geq tau, with tau>1.
k_star=min(find((stopCrit<=1)==1));
if ~isempty(k_star)&&k_star>1
    k_star=k_star-1;
end
noisyRMSE=norm(F_orig-F_data,'fro')/sqrt(m*n); %original RMSE error
noisySNR=20.*log10(norm(F_orig,'fro')/norm(F_orig-F_data,'fro'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maxIters=params(1); %time iterations in solving for wk
% dt=params(2);
% epsilon=params(3); %for regularizing TV
% lambda0=params(4); 
% q = params(5);
% alp0=params(6); %for tight/refined


%Plot the original, degraded, and recovered images
fullscreen=get(0,'ScreenSize');
font = 28; %fontsize setting
figure('Position',fullscreen)
image(F_orig); axis image; axis off; colormap(gray(256));
title('Original Image','FontSize',font)
% 
figure('Position',fullscreen)
image(F_data); axis image; axis off; colormap(gray(256));
title(['Noisy Image. RMSE=',num2str(noisyRMSE),', SNR=',num2str(noisySNR)],'FontSize',font)
if saveFlag==1
    figName=filePrefix+figPrefix+"noisy.png";
    saveas(gcf,figName)
    im=imread(char(figName));
    im_crop=crop_borders(im,[255],[0.01]);
    imwrite(uint8(im_crop(:,:,1)),char(figName));
    imwrite(uint8(F_data),char(filePrefix+"clean_"+figPrefix+"noisy.png"))
end

%Restored image
figure('Position',fullscreen)
image(xkArray(:,:,1,mink)); axis image; axis off; colormap(gray(256));
%subtitle="\lambda_k ="+lambda0+ "*"+q +"^k, \epsilon= "+epsilon+", \Delta t= "+ dt+", maxIters= "+maxIters;
titleStr="k="+num2str(mink-1)+", RMSE="+num2str(minVal)+", SNR="+num2str(snr(mink));
title(titleStr,'FontSize',font)
if saveFlag==1
    figName=filePrefix+figPrefix+"restored.png";
    saveas(gcf,figName)
    im=imread(char(figName));
    im_crop=crop_borders(im,[255],[0.01]);
    imwrite(uint8(im_crop(:,:,1)),char(figName));
    imwrite(uint8(xkArray(:,:,1,mink)),char(filePrefix+"clean_"+figPrefix+"restored.png"))
end

%k_star image
if (~isempty(k_star)) && (k_star~=mink)
    figure('Position',fullscreen)
    image(xkArray(:,:,1,k_star)); axis image; axis off; colormap(gray(256));
    %subtitle="\lambda_k ="+lambda0+ "*"+q +"^k, \epsilon= "+epsilon+", \Delta t= "+ dt+", maxIters= "+maxIters;
    titleStr="k^*="+num2str(k_star-1)+", RMSE="+num2str(rmse_final(k_star))+", SNR="+num2str(snr(k_star));
    title(titleStr,'FontSize',font)
    if saveFlag==1
        figName=filePrefix+figPrefix+"restored_kstar.png";
        saveas(gcf,figName)
        im=imread(char(figName));
        im_crop=crop_borders(im,[255],[0.01]);
        imwrite(uint8(im_crop(:,:,1)),char(figName));
        imwrite(uint8(xkArray(:,:,1,k_star)),char(filePrefix+"clean_"+figPrefix+"restored_kstar.png"))
    end
end


%Plot montage of multiscales near optimal
figure('Position',fullscreen)
numPlots =min(9,length(xkArray(1,1,1,:)));
montage(xkArray(:,:,:,end-numPlots+1:end),gray(256))

titleStr = "Multiscales k = "+num2str(numScales -numPlots)+" through k= "+num2str(numScales-1);
title(titleStr,'FontSize',font)

if saveFlag==1
    figName=filePrefix+figPrefix+"multiscales.png";
    saveas(gcf,figName)
    im=imread(char(figName));
    im_crop=crop_borders(im,[255],[0.01]);
    imwrite(uint8(im_crop(:,:,1)),char(figName));
end

%Plot montage of residuals near optimal
figure('Position',fullscreen)
numPlots =min(9,length(xkArray(1,1,1,:)));
residuals = xkArray;
residuals(:,:,1,:) = abs(F_orig-residuals(:,:,1,:))+125;
montage(residuals(:,:,:,end-numPlots+1:end),gray(256))

titleStr = "Residuals k = "+num2str(numScales -numPlots)+" through k= "+num2str(numScales-1);
title(titleStr,'FontSize',font)

if saveFlag==1
    figName=filePrefix+figPrefix+"residuals.png";
    saveas(gcf,figName)
    im=imread(char(figName));
    im_crop=crop_borders(im,[255],[0.01]);
    imwrite(uint8(im_crop(:,:,1)),char(figName));
end

%Plot the RMSE and SNR, as well as the stopping criterion plot
fig=figure('position',[100,100,1300,500]);
set(fig,'defaultAxesColorOrder',[[0,0,0];[0,0,0]]);
subplot(1,2,1)
yyaxis left
plot(0:numScales-1,rmse_final,'k-o','LineWidth',2)
xlabel('Multiscales: k','FontSize',22)
ylabel('RMSE','FontSize',22)
title(['RMSE & SNR vs multiscale-decompositions, kMin=',num2str(mink-1)],'FontSize',22)
yyaxis right
plot(0:numScales-1, snr,'r-.^','LineWidth',2)
ylabel('SNR','FontSize',22)
legend({'RMSE','SNR'},'FontSize',18)
subplot(1,2,2)
semilogy(0:numScales-1, stopCrit,'k-o','LineWidth',2)
hold on
semilogy(0:numScales-1,ones(numScales,1),'k--','LineWidth',2)
xlabel('Multiscales: k','FontSize',22)
if tightFlag(1)==1
    title(['(H(F_{data},Tx_k)+\alpha_kTV(log(x_k)))/H(F_{data},Tu), k^*=',num2str(k_star-1)],'FontSize',22)    
else
    title(['H(F_{data},Tx_k)/D(F_{data},u), k^*=',num2str(k_star-1)],'FontSize',22)
end

if saveFlag==1
    figName=filePrefix+figPrefix+"metrics.png";
    saveas(gcf,figName)
    im=imread(char(figName));
    im_crop=crop_borders(im,[255],[0.01]);
    imwrite(uint8(im_crop),char(figName));
end