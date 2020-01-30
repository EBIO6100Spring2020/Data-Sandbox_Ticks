# Data-sandbox
An initial space for data downloading and EDA.

Since most projects will share data sources, encounter similar issues, and need similar coding, use this repository to download data for Exploratory Data Analysis (EDA) for all the projects together. Once we have a clearer idea about data issues, we can split everything out into separate project repositories.

Use the `Issues` tab in this repository to post questions about data, make TO DO lists etc.

## Repo structure
There are many possible structures for a repository and [no real standards](https://www.explainxkcd.com/wiki/index.php/927:_Standards). Here is something to kick us off.

Subdirectories
* `data_raw` - put raw data here and never modify it
* `data_derived` - derived datasets (e.g. cleaned, sampled, reworked)
* `source` - any custom functions sourced by scripts
* `output` - temporary figures, tables, etc that you want to save
* `docs` - documentation and resources

Scripts (`.R`, `.Rmd`) can live at the top level. Name scripts by project keyword, function, and  perhaps include initials if you are working independently to start with (e.g. `ticks_eda_bam.R`). As a pipeline develops, include the sequence structure (e.g. `ticks01_download.R`, `ticks02_clean.R`, `ticks03_eda.R`,...).

## To set it up
Clone the repository to your computer using R studio as described in [Happy Git with R 12.3](https://happygitwithr.com/rstudio-git-github.html#clone-the-new-github-repository-to-your-computer-via-rstudio). See the [Git tutorials](https://github.com/EBIO6100Spring2020/Class-materials/tree/master/tutorials) if you need a refresher.

## Coding tips
### Refer to files in scripts using relative rather than absolute paths.
Use relative paths so that code will work from any location on any computer. Don't use absolute paths in scripts such as
```
C:/user/jane/janescoolstuff/experiment2/data_raw/neon_ants.csv
```
This will break the script on a different user's computer. Instead, use relative paths, such as
```
data_raw/neon_ants.csv
```
Anyone can then run the code without needing to modify the file paths. This is especially important when collaborating via a repository.

### Don't use `setwd()` in scripts.
This is not portable and will break the script on another person's computer. If you set up an RStudio project (as above), you will be in the correct working directory when you start RStudio by opening the project.

### Don't save your R workspace
[Start clean each time](https://r4ds.had.co.nz/workflow-projects.html#what-is-real). RStudio setup: In Tools > Global options > General, set "save workspace" to "Never" and uncheck everything except "Automatically notify me of updates to RStudio". This ensures that all your work derives from code and provides a test of the code each time you work on the script.
