# MHDM
This repository contains files to run MHDM decomposition on images affected by multiplicative gamma noise. It assumes Test Images is stored in the same directory. 
The scripts to run the tests are:
* testImagesAdditive.m
* testImagesAdditiveTight.m
* testImagesAdditiveRefined.m

They rely on the following functions
- Osher.m
- OsherTight.m
- OsherRefined.m
- metrics.m
- plotFigsOsher.m

Images are stored in the folders under the "Additive" main folder. The subfolders describe the contents. 

Files for the AA model with TV(log(u)) penalty term are complete for regular and tight versions. For the AA model, the test images have been slightly modified by adding 1 to all pixel values of the original images. Then, a fixed random gamma noise is applied (the random number generator is seeded with 10 for all runs). The scripts to run the tests are: 
* testImages.m
* testImagesTight.m

They rely on the following functions
- AAlog_blur.m
- AAlog_blur_tight.m
- metricsAA.m
- plotFigsAA.m

Output images are stored in the AA folder with appropriately labled subfolders. 
