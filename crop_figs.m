%script to crop figs so borders are trimmed.
folderPrefix=["AA/","additive/"];
folderNames=["barbara","cameraman","mandril","disc_square"]; 
cases=["regular","tight","refined"];
    
for pfx=folderPrefix
    for name=folderNames
        for cs=cases
            switch cs
                case "regular"
                    if pfx=="additive/"
                        folderPath=pfx+name+"_noise_additive/";   
                    else 
                        folderPath=pfx+name+"_noise/";
                    end
                case "tight"
                    folderPath=pfx+name+"_noise_tight/";
                case "refined"
                    folderPath=pfx+name+"_noise_refined/";
            end
            in_folder=char("./"+folderPath);
            out_folder=char("./crops/"+folderPath);
            p=0.01;
            store_crops(p,in_folder,out_folder);
        end
    end
end

function store_crops(p,in_folder,out_folder)
%performs crops on all images in in_folder, and stores them in out_folder
%with identical names to those in in_folder
%p=0.01; %percent to pad crop by.
%in_folder='./additive/cameraman_noise_additive/';
%out_folder="./crops/additive/cameraman_noise_additive/";

mkdir(char(out_folder))
files = loadImages(in_folder,'*.png');
for idx=1:length(files)
    figFile=files{idx};
    if ~contains(figFile,'clean')
        im=imread([in_folder figFile]);
        im_crop=crop_borders(im,[255],[p]);
        imwrite(uint8(im_crop(:,:,1)),char(out_folder+"crop_"+figFile));
    end
end

end
function Seq = loadImages(imgPath, imgType)
    %imgPath = 'path/to/images/folder/';
    %imgType = '*.png'; % change based on image type
    images  = dir([imgPath imgType]);
    N = length(images);
    % check images
    if( ~exist(imgPath, 'dir') || N<1 )
        disp('Directory not found or no matching images found.');
    end
    % preallocate cell
    Seq{N,1}=[];
    for idx = 1:N
        Seq{idx} = images(idx).name;
    end
end