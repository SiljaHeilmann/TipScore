# TipScore
**MATLAB** functions for image analysis of branched organ structures (such as the pancreas). See details [here](https://docs.google.com/document/d/143IZt6-4IdLZK5zwaFn9D988-nppJYnUX7sv5Kt3IIc/edit)

- Make sure that current Matlab directory contains the folder ``` sample_images ``` and all function .m files listed further below:

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

