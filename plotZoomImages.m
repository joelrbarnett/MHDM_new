%Create zoomed/detailed figures of images for publication
function plotZoomImages(F_orig, F_data, xkArray,zoomBox,filePrefix,figPrefix,saveFlag,tightFlag)

numScales = length(xkArray(1,1,1,:));
[~,rmse_final,stopCrit,~]= metrics(F_orig,F_data,squeeze(xkArray),numScales,tightFlag);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find min error and stopping criteria, original RMSE
%Min error and index of that error. Divide by (m*n) to get RMSE
[~,mink]=min(rmse_final);
%k_star = max_k D(F_data,Txk)^2/ D(F_data,Tu)^2 \geq tau, with tau>1.
k_star=min(find((stopCrit<=1)==1));
if ~isempty(k_star)&&k_star>1
    k_star=k_star-1;
end


%Plot the original, degraded, and recovered images at the specified zoom.
%Save actual zoomed regions as images
xl=zoomBox(1); xr=zoomBox(2); yl=zoomBox(3); yr=zoomBox(4);

figure()
image(F_orig(xl:xr,yl:yr)); axis image; axis off; colormap(gray(256));
title("Original Image")
if saveFlag==1
    figName=filePrefix+figPrefix+"zoom.png";
    imwrite(uint8(F_orig(xl:xr,yl:yr)),char(figName))
end
% 
figure()
image(F_data(xl:xr,yl:yr)); axis image; axis off; colormap(gray(256));
title('Noisy Image')
if saveFlag==1
    figName=filePrefix+figPrefix+"noisy_zoom.png";
    imwrite(uint8(F_data(xl:xr,yl:yr)),char(figName));
end

%Restored image
figure()
image(xkArray(xl:xr,yl:yr,1,mink)); axis image; axis off; colormap(gray(256));
titleStr="Restored: k="+num2str(mink);
title(titleStr)
if saveFlag==1
    figName=filePrefix+figPrefix+"restored_zoom.png";
    imwrite(uint8(xkArray(xl:xr,yl:yr,1,mink)),char(figName))
end

%k_star image
if (~isempty(k_star)) && (k_star~=mink)
    figure()
    image(xkArray(xl:xr,yl:yr,1,k_star)); axis image; axis off; colormap(gray(256));
    titleStr="Restored: k^*="+num2str(k_star);
    title(titleStr)
    if saveFlag==1
        figName=filePrefix+figPrefix+"restored_kstar_zoom.png";
        imwrite(uint8(xkArray(xl:xr,yl:yr,1,k_star)),char(figName))
    end
end


