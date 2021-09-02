# TipScore
**MATLAB** functions for image analysis of branched organ structures (such as the pancreas). See details: [STAR Protocols](https://www.cell.com/star-protocols/home)

- Download test data folders ``` raw_images_tif_files ``` and ``` outlines_mat_files ``` from [Mendeley Data](https://data.mendeley.com/drafts/nr9cyyk265). (574MB) 

- Make sure that current Matlab directory contains the folders ``` raw_images_tif_files ``` and ``` outlines_mat_files ``` and all function .m files listed further below:

- Run ``` 'demoScript_OneImage_ManualInput.m' ``` to see an example of a script that analyses one image with manual input from the user during the run.

- Run ``` 'demoScript_Batch_NoInput.m' ``` to see an example of a script that analyses all 7 sample images with no manual input from the user.

- Note that both scripts require **'Image Processing Toolbox'** and **'Statistics and Machine Learning Toolbox'** and 2018a version MATLAB or newer to run.

- Functions called by ``` demoScript_OneImage_ManualInput.m ``` and  ``` demoScript_Batch_NoInput.m ``` :

```
segmentDAPIimage.m
tipScoreIm.m
returnTableWithCellInt.m
returnTableWithCellPairInt.m
```


## Graphical abstract: 
<img width="450" alt="Screen Shot 2021-07-07 at 13 04 12" src="https://user-images.githubusercontent.com/11952601/124749116-1f302e00-df24-11eb-9c5c-f15845e7acd8.png">


**S. Heilmann, H. Semb, P. Nyeng.** *Quantifying spatial position in a branched structure in immunostained mouse tissue sections*. STAR Protocols, 2021.
