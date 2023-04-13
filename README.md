# qPCR-analysis from Cq values

Takes excel file output from qPCR machine, sample names layout and optional standard curve parameters. Plots Cq data, absolute concentrations and Tms in a neat format, ~ almost ready for publication.
- Code works with .xlsx file output by `Applied Biosystems Quantstudio-3`)
- If you need to analyze data from other instruments and if the current code fails, export to RDML format and use the script `linregpcr_processing.R`

_standard curve_ : Include Stdxx (xx = numbering) in the excel file name to process data and make a standard curve using data within the file

*If you're wondering that you can do all this with your expert excel skills, here's my rationale to convince you*: Coding this trivial analysis and plotting is very useful because I run tons of qPCR and I want all the plots to have the same formatting for easy comparision of data between runs without clicking around for every analysis. I wouldn't want to plot all these on excel! The script is especially handy when I run some experiment in a hurry for a meeting and I can get something to show in 2 mins by running the script 😂

## Instructions to run the tool
1. Clone the respository using git onto your computer. This brings in the code files and also mimics the folder structure (check if all the necessary folders are present)
2. Open the scripts in R by clicking on the `qPCR analysis.Rproj` file. _This ensures that R's current working directory is at the base of this project. 
	- If you haven't you need to install R, Rstudio, git; and install the packages listed with the library commands in `0-general_functions_main.R` script using `install.packages('..')`_
3. *Set metadata sheets*: Default option is to add your plate layout in the `excel files/Plate layout.xlsx` file. Make sure the `template_source` variable is set to *'excel'* in `0.5-user-inputs.R` file.
	- If you would prefer a cloud option, then make a googlesheet with sheets named `qPCR plate layouts` and `qPCR Std curves` following the template from [my sheet](https://docs.google.com/spreadsheets/d/1RffyflHCQ_GzlRHbeH3bAkiYo4zNlnFWx4FXo7xkUt8/edit#gid=0)
	- Copy the url of your googlesheet and replace the one in `sheeturls = list(plate_layouts_PK = ..` in the script `0-general_functions_main.R`
5. Open `0.5-user-inputs.R` script. Change the `flnm` and `title_name` variables. and other options by following the comments in the file
6. You can try running the template file `q33_RNA stability-3_20-6-22` within the `excel files` folder by copying the template data for q33 and Std30 standard data from [my sheet](https://docs.google.com/spreadsheets/d/1RffyflHCQ_GzlRHbeH3bAkiYo4zNlnFWx4FXo7xkUt8/edit#gid=0)
7. Open `analysis.R` script and source it to analyze the data

**Sample naming scheme**
Sample names are written in a google sheet (contact me for access to see the template) in this particular format
`targetname-samplecategory_templatename.biologicalreplicatenumber`
- targetname = name of the amplicon, _ex: 16s, U64, spliced_
- samplecategory could stand for test, positive, negative (controls) or NoRT control etc.
- templatename = name of the DNA or RNA sample loaded as the template
- biologicalreplicatenumber = 1,2,3 ..
Examples: old16sU64-Test_295.1,	old16sU64-NoRT_89.1

![image](https://user-images.githubusercontent.com/14856479/113488074-6cae8200-9481-11eb-9d82-e97033b72e2e.png)

The output image looks like this
- primerpairname becomes facets or panels (boxes on the top); 
- template on the x axis (named `assay_variable`) and 
- samplecategory as colour (named as `Sample_name`)
<img src = 'https://user-images.githubusercontent.com/14856479/113488826-1859d100-9486-11eb-8384-1ad17afea737.png' width = "600">



**Git organization**

1. *main* branch holds the latest version of the code for general purpose plotting for any qPCR output file with sample names in the above way
2. For customized plots for individual experiments - there will be an ad-hoc branch with specific code files that modify these base plots. Or for publications there will be a different branch named after the experiment (ex: S024) with the highly customized plotting formats and extra items/tables output for supplmentary data (this will make it easy to provide the code along with the publication)
3. There are different main branches for each of the major tasks the experiments fall under a. Time_series_master (ex: S019) b. master (ex: S024) (for generic qPCR) C. Dose_response (ex: S03) : If branches with these names don't exist yet, you can start with the branches named after the experiments


# qPCR analysis using LinRegPCR

Linregpcr tool analyzes raw qPCR curves and calculates the ~initial fluorescence which can be converted into initial concentration using calibration standards (6 replicates of the same concentration known standard). I have made a wrapper to run their python version of LinRegPCR through R, but the portability and documentation is in progress.

If you want to try it out, the script is called `linregpcr_processing.R`

*[Source](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04306-1): Untergasser, Andreas, et al. "Web-based LinRegPCR: Application for the visualization and analysis of (RT)-qPCR amplification and melting data." _BMC bioinformatics_ 22.1 (2021): 1-18.)*.

# Copyleft : GPL-3.0-or-later license
```
wrappers for automated processing and plotting of bacterial flow cytometry data 
Copyright (C) 2023  Prashant Kalvapalle

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
```
