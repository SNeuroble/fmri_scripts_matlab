# Stephanie's fMRI Analysis Repository - Matlab

The contents are standalone scripts for fMRI analysis and a few other general tasks in Matlab.

## Getting Started

### Prerequisites

The only prerequisite is a working version of Matlab (haven't been tested yet across versions, but everything is expected to work with Matlab 2018b or greater).

### Overview

A few of these materials have been developed to conduct the core analyses for a project and are described below.

#### 1. Generalizability Theory ICC toolbox

*Location:* *reliability > ICC > ICC\_toolbox*

*Purpose:* Calculate ICC coefficients including multiple facets (sources of error).

*Usage:*

1. Add all toolbox files to your path.

2. Load data and factor table. This can be accomplished using the *load\_reliability\_data* script, which automatically creates a factor table from the filenames. This can be changed for your data or data can be loaded manually (try "help load\_reliability\_data" for more info on the format). 

3. Run reliability analysis. Use *run\_reliability*. Try "help run\_reliability" for usage.

*References:* 

Noble, S., Spann, M. N., Tokoglu, F., Shen, X., Constable, R. T., & Scheinost, D. (2017). Influences on the test–retest reliability of functional connectivity MRI and its relationship with behavioral utility. Cerebral Cortex, 27(11), 5415-5429.

Noble, S., Scheinost, D., Finn, E. S., Shen, X., Papademetris, X., McEwen, S. C., ... & Mirzakhanian, H. (2017). Multisite reliability of MR-based functional connectivity. Neuroimage, 146, 959-970.

#### 2. Summarizing true positives in cluster extent-based inference

*Location:* *clf*

*Purpose:* Summarization and light analysis of true positive maps previously calculated elsewhere. Bash scripts for estimating true positives from task data can be found in the *power\_cluster\_failure repository*.

(Reference forthcoming) 



## Contact

Please let me know if you have any questions about the materials provided here! Although I try to make all my code publicly available, some of the materials in this repository have been developed for more private use and may therefore require further explanation.



