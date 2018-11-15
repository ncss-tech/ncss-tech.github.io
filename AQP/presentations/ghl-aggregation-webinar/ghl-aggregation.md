<!-- 
ideas:
http://landeco2point0.wordpress.com/2013/11/07/r-studio-and-presentations-and-git-oh-my/
http://rstudio-pubs-static.s3.amazonaws.com/27777_55697c3a476640caa0ad2099fe914ae5.html
https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/aqp/gen-hz-assignment.html?root=aqp
-->





Soil Data Aggregation using R
========================================================
transition: none
width: 1024
height: 800
css: custom.css

<br><br><br><br><br><br><br><br><br>
<span style="color: white; font-size:75%;">This document is based on:</span>
<ul style="color: white; font-size:75%;">
<li>`aqp` version 1.8</li>
<li>`soilDB` version 1.5-2</li>
<li>`sharpshootR` version 0.7-2</li>
</ul>

What is R?
================
- an analysis platform: calculator, statistics program, GIS, etc...
- R syntax: vocabulary to explore, summarize, and model data
![alt text](static-figures/rstudio.png)


What is R?
================
- an integrator: analysis + GIS + database connectivity
- ODBC and GDAL link R to nearly all possible formats/interfaces
![alt text](static-figures/triangle.png)


What can R do?
========================================================
- Base R (*functionality is extended through packages*)
  - basic summaries of quantitative or qualitative data
  - data exploration via graphics
  - GIS data processing and analysis

- Integrated Workflow (*self-documenting, repeatable*)
  - report writing: code + summaries in a single document
  - digital soil mapping: modeling properties and classes
  - evaluation of OSD, KSSL, NASIS, ... for SDJR or related

- Soil Science Specific Packages
  - visualization, aggregation, classification <span class="link-to-details">&#8594;&nbsp;[aqp package](http://aqp.r-forge.r-project.org/aqp-html-manual/index.html)</span>
  - access to commonly used soil databases <span class="link-to-details">&#8594;&nbsp;[soilDB package](http://aqp.r-forge.r-project.org/soilDB-html-manual/index.html)</span>
  - misc. soil survey specific functions <span class="link-to-details">&#8594;&nbsp;[sharpshootR package](http://aqp.r-forge.r-project.org/shapshootR-html-manual/index.html)</span>
 

What can R do?
========================================================

<br><br>
![alt text](static-figures/loafercreek-ml-hz.png)
<br><br>
<center>**estimate a most-likely horizonation, details to follow**</center>

What can R do?
========================================================

<br><br>
![alt text](static-figures/loafercreek-ric.png)
<br><br>
<center>**estimate reasonable RIC, details to follow**</center>

Setup your local machine
========================================================

- Be sure to follow embedded links to full documentation
- Check out the [AQP project website](http://aqp.r-forge.r-project.org/)
- Browse some of the [various documentation](http://cran.r-project.org/manuals.html) on R

- The first time you open RStudio, **install** required packages:

```r
install.packages('aqp', dep=TRUE) 
install.packages("RODBC", dep=TRUE)
install.packages('soilDB', dep=TRUE)
```

- Every subsequent time, **load** required packages:

```r
library(aqp)
library(soilDB)
library(latticeExtra)
# ... other libraries you might be using
```


Connect to NASIS via ODBC
========================================================

- connection allows R to pull data from a NASIS *selected set* via the soilDB package
- connection sets up a small file that stores the NASIS login and password
- this is a one-time operation <span class="link-to-details">&#8594;&nbsp;[ODBC setup instructions](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/soilDB/setup_local_nasis.html?root=aqp)</span>

- getting data from NASIS is as simple as:

```r
# get pedons from the selected set
pedons <- fetchNASIS()

# plot the first 15
plot(pedons[1:15, ])

# get component data from the selected set
components <- get_component_data_from_NASIS_db()
```
- be sure to read the documentation <span class="link-to-details">&#8594;&nbsp;[`fetchNASIS()`](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/soilDB/soilDB-Intro.html?root=aqp)</span>


Aggregation of Pedon Data: Common Problems
========================================================
![alt text](static-figures/genhz-sketch.png)

- description styles
 - the age-old splitters vs. lumpers
 - use of transitional horizons - I say tomato, you say tomato!
  
- legacy horizon nomenclature
 - old school meets new school and they ain't speakin' the same language.....

- soils described to varying depths - it's the pits!


Generalized Horizon Labels: Micro-correlation
========================================================
![alt text](static-figures/genhz-sketch.png)

- generalized horizon labels = GHL
- determine the core concept 
- which way do transition horizons go?
 - GHL process allows for flexibility
- horizon redundancy - how many Bt horizons are necessary?
- outlier horizons - what happens to them?


Summary of GHL workflow
========================================================

<center>**Most of these steps are implemented in the reports presented last**</center>

1. Select a set of GHL that best represents a group of pedons to be aggregated.

2. Assign GHL to each horizon using whatever information is available for grouping horizons.

3. Evaluate GHL assignments and manually refine as needed.

4. Keep track of final GHL assignments in NASIS or text file.

5. Estimate a most likely horizonation e.g., top and bottom depths for each generalized horizon label.

6. Compute range in characteristics, aka low-rv-high values, for clay, sand, pH, etc. over GHL. (next document in this series)

[&#8594;&nbsp;extended tutorial](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/aqp/gen-hz-assignment.html?root=aqp)



Assumptions
========================================================
- a collection of pedons for a component should represent an aggregated component
- a GHL process of some kind is required to efficiently summarize and aggregate horizon data from these collections
- low-rv-high values must be tied to some metric: 
 - 5th-50th-95th percentile is one option
- linked pedons represent examples of the aggregate data, but aggregate data should not be limited to a very small set of linked pedons

- **IMPORTANT** these reports and data summaries are a useful *starting point* for the aggregation process




Summarizing Data with Quantiles (Percentiles)
========================================================

- suggested basis for RIC: 
  - 5th, 50th, and 95th percentiles are low-rv-high values
  - no assumptions of distribution, simple interpretation

<img src="ghl-aggregation-figure/quantiles-normal-dist-1.png" title="plot of chunk quantiles-normal-dist" alt="plot of chunk quantiles-normal-dist" style="display: block; margin: auto;" />


Summarizing Data with Quantiles (Percentiles)
========================================================

- suggested basis for RIC: 
  - 5th, 50th, and 95th percentiles are low-rv-high values
  - no assumptions of distribution, simple interpretation

<img src="ghl-aggregation-figure/quantiles-exp-dist-1.png" title="plot of chunk quantiles-exp-dist" alt="plot of chunk quantiles-exp-dist" style="display: block; margin: auto;" />


A Sample Dataset
=======================================================

- 15 pedons correlated to the [Loafercreek](https://soilseries.sc.egov.usda.gov/OSD_Docs/L/LOAFERCREEK.html) soil series
- Fine-loamy, mixed, superactive, thermic Ultic Haploxeralfs
- Common soil formed on meta-volcanic rocks of the Sierra Nevada Foothills, MLRA 18
- Included in the `soilDB` package for testing purposes

<img src="ghl-aggregation-figure/load-data-1.png" title="plot of chunk load-data" alt="plot of chunk load-data" style="display: block; margin: auto;" />

Consult the OSD for Ideas
========================================================
- Look up the series RIC if available

![alt text](static-figures/RIC.png)

- Horizons from the OSD may be a good starting point for GHL template



Tabulate Horizon Designations
========================================================

- sort by frequency

|  A| Bt1| Bt2| Bt3| Oi|  R| Cr| Crt| BA| BCt| Bt4| Bw| CBt| 2BCt| 2Bt3| 2CB| 2Cr| 2R|
|--:|---:|---:|---:|--:|--:|--:|---:|--:|---:|---:|--:|---:|----:|----:|---:|---:|--:|
| 15|  15|  15|  10|  9|  9|  8|   5|  3|   3|   2|  2|   2|    1|    1|   1|   1|  1|

- sort alphabetically

| 2BCt| 2Bt3| 2CB| 2Cr| 2R|  A| BA| BCt| Bt1| Bt2| Bt3| Bt4| Bw| CBt| Cr| Crt| Oi|  R|
|----:|----:|---:|---:|--:|--:|--:|---:|---:|---:|---:|---:|--:|---:|--:|---:|--:|--:|
|    1|    1|   1|   1|  1| 15|  3|   3|  15|  15|  10|   2|  2|   2|  8|   5|  9|  9|

- plot ranges in horizon depths
<img src="ghl-aggregation-figure/horizonation-mid-point-1.png" title="plot of chunk horizonation-mid-point" alt="plot of chunk horizonation-mid-point" style="display: block; margin: auto;" />


Summarize Available Soil Properties
========================================================


<img src="ghl-aggregation-figure/univariate-eval-clay-1.png" title="plot of chunk univariate-eval-clay" alt="plot of chunk univariate-eval-clay" style="display: block; margin: auto;" />

<img src="ghl-aggregation-figure/univariate-eval-rf-1.png" title="plot of chunk univariate-eval-rf" alt="plot of chunk univariate-eval-rf" style="display: block; margin: auto;" />


Determination of a GHL Template and Rules
========================================================
- pattern matching via [regular expression](http://www.regexr.com/) (REGEX)
 - this is where most micro-correlation decisions are defined

- GHL and rules for our sample dataset:
  - **A**: `^A$|Ad|Ap`
  - **Bt1**: `Bt1$`
  - **Bt2**: `^Bt2$`
  - **Bt3**: `^Bt3|^Bt4|CBt$|BCt$|2Bt|2CB$|^C$`
  - **Cr**: `Cr`
  - **R**: `R`

- special characters in REGEX rules: 
 - `|` = "or"
 - `^` = anchor to left-side
 - `$` = anchor to right-side






Assignment of GHL
========================================================

- cross-tabulation of original names (columns) vs. GHL (rows)
- be sure to check the `not-used` row


|         | 2BCt| 2Bt3| 2CB| 2Cr| 2R|  A| BA| BCt| Bt1| Bt2| Bt3| Bt4| Bw| CBt| Cr| Crt| Oi|  R| Sum|
|:--------|----:|----:|---:|---:|--:|--:|--:|---:|---:|---:|---:|---:|--:|---:|--:|---:|--:|--:|---:|
|A        |    0|    0|   0|   0|  0| 15|  0|   0|   0|   0|   0|   0|  0|   0|  0|   0|  0|  0|  15|
|Bt1      |    0|    0|   0|   0|  0|  0|  0|   0|  15|   0|   0|   0|  0|   0|  0|   0|  0|  0|  15|
|Bt2      |    0|    0|   0|   0|  0|  0|  0|   0|   0|  15|   0|   0|  0|   0|  0|   0|  0|  0|  15|
|Bt3      |    1|    1|   1|   0|  0|  0|  0|   3|   0|   0|  10|   2|  0|   2|  0|   0|  0|  0|  20|
|Cr       |    0|    0|   0|   1|  0|  0|  0|   0|   0|   0|   0|   0|  0|   0|  8|   5|  0|  0|  14|
|R        |    0|    0|   0|   0|  1|  0|  0|   0|   0|   0|   0|   0|  0|   0|  0|   0|  0|  9|  10|
|not-used |    0|    0|   0|   0|  0|  0|  3|   0|   0|   0|   0|   0|  2|   0|  0|   0|  9|  0|  14|
|Sum      |    1|    1|   1|   1|  1| 15|  3|   3|  15|  15|  10|   2|  2|   2|  8|   5|  9|  9| 103|



Evaluation of GHL
========================================================

- plot profiles with horizons colored according to GHL
- does it make sense?
- what about horizons in the `not-used` group?

<img src="ghl-aggregation-figure/plot-ghl-1-1.png" title="plot of chunk plot-ghl-1" alt="plot of chunk plot-ghl-1" style="display: block; margin: auto;" />


Evaluation of GHL
========================================================
- plot range in GHL depths, check for overlap

<img src="ghl-aggregation-figure/eval-ghl-2-1.png" title="plot of chunk eval-ghl-2" alt="plot of chunk eval-ghl-2" style="display: block; margin: auto;" />


Evaluation of GHL
========================================================
- multivariate summary of depth, clay content, and RF volume



<img src="ghl-aggregation-figure/mds-plot-1.png" title="plot of chunk mds-plot" alt="plot of chunk mds-plot" style="display: block; margin: auto;" />

Importing GHL to NASIS
========================================================
- generalized horizon labels (so far) are stored in R session
- Use an R script to create a text file, called `horizon_agg.txt`, of horizon ID's and corresponding GHL assignments
- run the NASIS calculation to update the comp layer id field with the GHL for each pedon horizon
- Calculations/Validations --> Pedon Horizon 
  - 'Update horizon group aggregations using a text file'
  - Calculation will look for the `horizon_agg.txt file` in: `C:/data/horizon_agg.txt`
  - Calculation will load the file and update the comp layer id field
- manual adjustments in NASIS: outliers, special cases, etc.

[&#8594;&nbsp;full instructions here](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/aqp/gen-hz-assignment.html?root=aqp)




Demonstration of Several Report Styles
========================================================

- Stephen's Examples: 
  - [Lecyr Pedon Report](https://github.com/sroecker01/soil-pit/blob/master/examples/lecyr_pedon_report.md)
  - [Genesee Lab Data Report](https://github.com/sroecker01/soil-pit/blob/master/examples/genesee_lab_report.md)
  - [Cincinanti Map Unit Report](https://github.com/sroecker01/soil-pit/blob/master/examples/cincinnati_mapunit_report.md)

- Dylan's Examples:
  - [Loafercreek Pedon Report](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/example-reports/loafercreek.html?root=aqp)
  - [Dunstone Pedon Report](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/example-reports/dunstone.html?root=aqp)
  - [Amador Pedon Report](https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/example-reports/amador.html?root=aqp)


Thank You
========================================================
**Questions, comments, ideas**
<br><br>

**R resources for pedologists**
- <span class="link-to-details">&#8594;&nbsp;[NCSS Job-Aids](http://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/edu/ncss/?cid=nrcs142p2_054322)</span>
- <span class="link-to-details">&#8594;&nbsp;[aqp tutorials](http://aqp.r-forge.r-project.org/)</span>
- <span class="link-to-details">&#8594;&nbsp;[Statistical data analysis for pedologists](http://www2.gru.wvu.edu/~tdavello/files/stats/table_of_contents.html)
- <span class="link-to-details">&#8594;&nbsp;[Dylan Beaudette's blog](http://casoilresource.lawr.ucdavis.edu/blog/)
- <span class="link-to-details">&#8594;&nbsp;[soil-pit Github repository](https://github.com/sroecker01/soil-pit)
<br><br>

**Additional AQP Contributors:**
- Pierre Roudier (Landcare Research)


