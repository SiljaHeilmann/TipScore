# TipScore
**MATLAB** functions for image analysis of branched organ structures (such as the pancreas). See details [Cell Press STAR Protocols](https://docs.google.com/document/d/143IZt6-4IdLZK5zwaFn9D988-nppJYnUX7sv5Kt3IIc/edit)

- Download test data folders ``` raw_images_lsm_files ``` and ``` outlines_mat_files ``` from [Mendeley Data](https://data.mendeley.com/drafts/nr9cyyk265). (580MB) 

- Make sure that current Matlab directory contains the folders ``` raw_images_lsm_files ``` and ``` outlines_mat_files ``` and all function .m files listed further below:

- Run ``` 'demoScript_OneImage_ManualInput.m' ``` to see an example of a script that analyses one image with manual input from the user during the run.

- Run ``` 'demoScript_Batch_NoInput.m' ``` to see an example of a script that analyses all 8 sample images with no manual input from the user.

- Note that both scripts require **'Image Processing Toolbox'** and **'Statistics and Machine Learning Toolbox'** and 2018a version MATLAB or newer to run.

- Functions called by ``` demoScript_OneImage_ManualInput.m ``` and  ``` demoScript_Batch_NoInput.m ``` :

```
segmentDAPIimage.m
tipScoreIm.m
returnTableWithCellInt.m
returnTableWithCellPairInt.m
```


#Graphical abstract: 

<img width="500" alt="Screen Shot 2021-04-13 at 13 40 59" src="https://user-images.githubusercontent.com/11952601/114546969-164dfa00-9c5e-11eb-9051-7f91bb9edbd6.png">

S. Heilmann, H. Semb, P. Nyeng. Quantifying spatial position in a branched structure in immunostained mouse tissue sections. STAR Protocols, Volume X, Issue X, 2021.
