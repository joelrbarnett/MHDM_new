%Script to create and save noisy images from the test image folder
clear all

%for saving
filePrefix="Test Images/";
figPrefix=["barbara_noise_02.png","cameraman_noise_02.png",...
    "pollen_noise_02.png","mandril_noise_02.png","circles_noise_02.png",...
    "geometry_noise_02.png","disc_square_02.png"];
figPrefix04=["barbara_noise_04.png","cameraman_noise_04.png",...
    "pollen_noise_04.png","mandril_noise_04.png","circles_noise_04.png",...
    "geometry_noise_04.png","disc_square_04.png"];%for standard deviation 0.4

fileNames=["barbara","cameraman","pollen","mandril","circles","geometry","disc_square"]; 
images=["barbara.png","cameraman.tif","pollen.tif","mandril_gray.tif","circles.tif","geometry.tif","disc_square.png"];
imagesPNG=["barbara.png","cameraman.png","pollen.png","mandril.png","circles.png","geometry.png","disc_square.png"];

folder_path="Test Images/";
save_path ="Test_Images_plus1/";


%Form noisy image: 
%%% Gamma noise %%%   
rng(10);
a=25;%6.25; %gamma noise with mean 1, standard deviation 0.2=1/sqrt(a). 
GamNoise=gamrnd(a,1/a,size(F_orig));
for k= 1:length(images)
    %choose images
    F_orig=imread(char(folder_path+images(k))); 
    F_orig=F_orig+1;
    imwrite(F_orig,char(save_path+imagesPNG(k)));%to convert tifs to png
    
    %form noisy images
    F_orig=double(F_orig);
    F_data=F_orig.*GamNoise; %multiply noise into image
    F_data=uint8(F_data);
    imwrite(F_data,char(save_path+figPrefix(k)));
end



