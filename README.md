# pMFMC
parametric multi-fidelity Monte Carlo estimation methods

## Supplementary Code Submission

This supplementary package includes R scripts, generated figures, and a technical derivation document that support the results in the manuscript. Each script corresponds to a specific experiment or case study described in the paper. All files are anonymized for peer review.

## Prerequisites

### R Version
This code has been tested with R version 4.0.0 or higher.

### Required R Packages
Before running the scripts, install the following R packages:

```r
install.packages(c(
  "tidyverse",    # Data manipulation and visualization
  "ggplot2",      # Plotting (included in tidyverse)
  "evd",          # Extreme value distributions
  "MASS",         # Statistical functions
  "mvtnorm",      # Multivariate normal distribution
  "qqplotr",      # Enhanced Q-Q plots
  "stats"         # Statistical functions (base R)
))
```

### Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/mjkim1001/pMFMC.git
   cd pMFMC
   ```

2. Install required packages using the command above

3. Run the R scripts in the `scripts/` directory

## Folder Structure
```
├── data/                           # Real data used for ship application
├── scripts/                        # R scripts for each analysis
│   ├── BernoulliAnalysis.R         # Bernoulli model
│   ├── Gumbel_muOnly.R             # Gumbel: estimation of mu only
│   ├── Gumbel_jointEstimation.R    # Gumbel: estimation of mu and sigma jointly
│   └── Ship_application.R          # Application to ship data
├── images/                         # Generated PNG figures
└── derivations/                    # Technical derivation document
    └── Gumbel_supplementary.pdf    # Derivations to help understand the code for Gumbel estimation
```

## Usage
Each script can be run independently. Navigate to the `scripts/` directory and run:

```r
# For Bernoulli analysis
source("BernoulliAnalysis.R")

# For Gumbel estimation (mu only)
source("Gumbel_muOnly.R")

# For joint Gumbel estimation
source("Gumbel_jointEstimation.R")

# For ship application
source("Ship_application.R")
```

## Notes
- All code is self-contained and does not rely on any external data files beyond what's included in the `data/` directory.
- Each R script automatically saves the corresponding figures in `.png` format to the `images/` folder.
- All plots are reproducible by running the scripts in the `scripts/` folder.
- Scripts should be run from the project root directory to ensure correct relative paths.
