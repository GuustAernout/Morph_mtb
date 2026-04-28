# Morph_mtb
THis repo ontains the development of the morph app which is used to predict protein family members in mycobacterium tuberculosis.
This shiny app is based on the morph algorithm. It was originally deigned for udse in plant. This app uses this algorithm to predict previously unknown protein pathway members in *Mycobacterium tuberculosis*. 


  

In this repository, two folders where attached named Data and MorphMtb.

A separate R script called allMtbPathways.R is also present.

## Data folder
In the data folder all pre-processed and app implemented data files can be found. The 6 included datasets are:
- ESX.txt
- clark.txt
- drug.txt
- inaki.txt
- primary.txt
- timecourse.txt
The results from the K-means and SOM clustering are also found in this folder.

## MorphMtb folder
4 scripts can be found in this folder which are crucial for the shiny app:  
1. **MorphMtb.R**  
    This R script contains the shiny interface code. This file is used to start the application.
2. **functions.R**  
    This R script calculates and stores all scores used for the analysis.
3. **rentrez.R**  
    In this R script, a connection with NCBI is made to retrieve information about the candidate genes.
4. **styles.css**  
    The css file determines look of the shiny GUI.

All of the data needed for the calculations is also present in this folder. 
The files in this folder are read and used in the calculations. Some files will be written in this folder when chosen to implement own expression data.
      

## Using the Shiny web application

To make use of the Shiny web application, R and RStudio need to be installed first. 

### Installing R on Windows

Go to the website using [this link](https://cran.r-project.org/bin/windows/base/) to install version 4.5.3 (the latest version at the moment) for Windows.

Click on the "Download R-4.5.3 for Windows" link. A download should be started.  
1. Open the downloaded file
2. Choose English as installation language
3. Start the setup
4. Do not customize the startup options
5. Finish setup

### Installing RStudio on Windows
Go to the website using [this link](https://posit.co/download/rstudio-desktop/) to install version 2024.04.1+748 (the latest version at the moment) for Windows.  

Click on the "Download RStudio desktop for Windows" button. A download should be started.  
1. Open the downloaded file
2. Click on install
3. Finish setup  
4. You can now open the RStudio application

### The MORPH-Mtb Shiny web application
To make use of the MORPH-Mtb Shiny web application, you first have to open the MorphMTB.R script with RStudio (Right-click on the file and choose open with RStudio).   
You will see a button with a green arrow and "Run App". Click on this button and the Shiny application will be started. You can choose to optionally open this in a web browser by clicking the "Open in Browser" button. 

The web application consists out of 3 pages:
1. **Gene centric query**  
    You can use this page if you only want to use the implemented 6 datasets.  
2. **Implement own expression data**    
    You can use this page if you want to upload your own expression data and use this together with the 6 datasets for the analysis.  
3. **About**  
    This is an informative page about the MORPH-Mtb algorithm. 


#### Gene centric query
ON this page queries can be performed solely on the 6 implemented daatsets. A target pathway is needed as input for the algorithm. You can input these via the text box or by uploading a text file. IMPORTANT: gene ids need to be enter separated.

Further, you can choose the number of random pathways generated, having the same length as the input pathway and the number of candidate pathways to display.  

After clicking the start button, your input pathway will be given in a table, telling you if there are gene IDs in your pathway that are not of *Mycobacterium tuberculosis*.  
In the next tab, Result random pathway, the AUSR-score of the random pathways is given and the average AUSR-score.  
In the last tab, Result input pathway, the ranked candidate genes are given together with some extra information:  
- The AUSR-score of your input pathway
- The dataset and the clustering solution that gave the highest score
- A table with the AUSR-scores from the datasets with their clustering solutions
- A table with the ranked candidate genes. Each gene gets a score and a description from NCBI.

This list of candidate genes can be downloaded in a .txt file by clicking the download link underneath the start button.  
The application can be restarted by clicking the "New analysis" button.  

Because the calculations can take some time, a spinning flower is shown while the application is working.

#### Implement own expression data
THis page cna be used if yo uwant to use yourown expression data. The data should contain the sample names as column names and the genes as row names. The expression counts need to be enter separated. You can take a look at the files of the already integrated files for an example.


The elbows of the plotted k-means clustering and SOM clustering will have to be defined. Further, you can choose the number of random pathways generated, having the same length as the input pathway and the number of candidate pathways to display.

After clicking the start button, the expression data will be given in a table and a processed file will automatically be written to the folder you're working in. This file will be used further down in the analysis.  
In the next tab, Filtered expression data, the percentage of genes kept after filtering will be given, so that you can decide if you want to work further with this data.  
In the third tab, Clustering, two plots (k-means and SOM clustering) will be visualed. Based on these plots, you can change the elbow and two clustering files will automatically be written in the folder you're working in.  
The next three tabs are the same as these of the first page and will give you the same output. The calculations will be done using the new expression data.  

The application can, again, be restarted by clicking the "New analysis" button.    
The spinning flowers indicate if the algorithm is running.

## allMtbPathways.R script

This R script was made for advanced bioinformaticians and people who know how to code. You can run this script with all of the Mtb pathways and get the candidate genes and AUSR score of each of them as output.  

The input of this script is an Excel file with on each sheet a different pathway of Mtb. This data will be written in new tab-separated files, which are then used to calculate the AUSR-score and the top 6 candidate genes. This will again be written in a new file so that the data is easily accessible.  
This way, all of the pathways of Mtb can be analyzed at the same time.  

You will have to adjust some code to get the script to run. First of all, the path to the folder where your excel file is stored on your laptop and then what the file is called. At the moment, only the first 6 candidate genes will be generated, but this can be altered.   

Further down the script, the top 3 AUSR-scored pathways are used to get the candidate genes and their description from. Out of these candidate genes, the ones with an unknown function are filtered and put into a new data frame. Research on these genes can be done and new information about Mtb can be found.  
Again, the amount of given candidate genes can be adjusted.

