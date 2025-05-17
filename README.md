# GAMMAmotif
**GAMMAmotif** is an R package for discovering gapped motifs using the GAMMA model. It provides a core function `GAMMA_fun()` and a built-in dataset `sample_data` to help users quickly get started with motif discovery.

## Installation
You can install the released version of demotif from GitHub with:
```
devtools::install_github("RanLIUaca/GAMMAmotif")
library(GAMMAmotif) 
```

## Prior settings
We set the non-informative prior distributions as the default settings for users: $\boldsymbol{\alpha_0}$ and $\boldsymbol{\alpha_jk}$ are two vectors whose entries are all one. $\boldsymbol{\pi_{0,a_i}}$ is set to be equal. $p_0$ is set to be 0.5. $\beta_{kj}$ and $\nu_{kj}$ are set to be 1.

## Input

For the `GAMMA_FUN` function, the input is:

- **`Data`**: A data frame containing two columns:
  - `seq`: A character vector of biological sequences.
  - `W_obs`: A numeric vector of observed motif labels. Use `NA` for unknown labels.

- **`dict`**: A character vector representing the alphabet (e.g., `["A", "C", "G", "T"]` for DNA).

- **`num_motif`**: *Integer*. The number of motifs to infer.

- **`Motif_length`**: *Integer vector*. The length of each motif (should be of length `num_motif`).

- **`sum_N`**: *Integer*. Total number of MCMC iterations to perform.

- **`burnin`**: *Integer*. Number of burn-in iterations before Metropolis-Hastings shift proposals are activated.

- **`mh_step_w`**: *Integer*. Frequency (in iterations) at which MH shift proposals are triggered.

---

##  Output

The `GAMMA_FUN()` function returns a **list** with the following components:

- **`Theta_0`**:  
  MAP (maximum a posteriori) estimate of the **background nucleotide distribution**, representing non-motif regions in the sequences.

- **`Theta`**:  
  A list of MAP-estimated **motif PWMs (position weight matrices)** for each inferred motif. Each PWM describes the base probabilities at each position within the motif.

- **`W`**:  
  MAP estimate of **motif assignments** for each sequence, indicating which motif (if any) is present.

- **`A`**:  
  A list of MAP-estimated **motif position matrices**, where each matrix specifies the inferred start positions of motifs across sequences.

- **`plot`**:  
  A list of `ggplot2` objects, including:
  - **Log-likelihood trace plot**: Shows convergence and stability of the MCMC algorithm.
  - **Sequence logo plots**: Visual representations of the discovered motif PWMs using `ggseqlogo`.



## Example
This is a toy example for decomposing hierarchical motifs. We shall use the data contained in the package.
``` R
# Load the GAMMAmotif package (contains the main inference function)
library(GAMMAmotif)

# Load required dependencies
library(MCMCpack)  
library(stringr)    
library(MASS)       
library(ggseqlogo)  
library(ggplot2)    

# Load the built-in example dataset
# The dataset contains 3 motifs, each of length 4, embedded in the sequences
data(sample_data)

# Preview the first few rows of the data
head(sample_data)

# Number of sequences
total_n = length(sample_data[,1])


# Define the alphabet (dictionary) used in the sequences
# Here we use the 20 amino acids for protein sequences
dict = 'ACDEFGHIKLMNPQRSTVWY'

# Split the string into individual characters
tmp_index = 1:str_length(dict)
dict = str_sub(dict, tmp_index, tmp_index)

# Create a mapping between characters and numeric indices
digit_dict = as.character(1:length(dict))
names(digit_dict) = dict
names(dict) = digit_dict

# Get dictionary length
len_dict = length(dict)

# -------------------------------
# Set inference parameters
# -------------------------------

num_motif = 3                  # Number of motifs to discover
Motif_length = rep(4, 3)       # Each motif is 4 characters long
sum_N = 20                     # Total number of MCMC iterations
burnin = 10                    # Burn-in period (used to discard early samples)
mh_step_w = 10                 # Frequency of Metropolis-Hastings "shift" proposals

# -------------------------------
# Run Bayesian motif inference
# -------------------------------

# GAMMA_FUN performs MCMC-based inference for motif PWM, positions (A),
# motif type assignments (W), and optional gap parameters (lambda)
Res = GAMMA_FUN(
  Data = sample_data,          # Input data: a data frame with 'seq' and 'W_obs'
  dict = dict,                 # Alphabet used in sequences
  num_motif = num_motif,       # Number of motifs to infer
  Motif_length = Motif_length, # Lengths of each motif
  sum_N = sum_N,               # Total MCMC iterations
  burnin = burnin,             # Number of burn-in iterations
  mh_step_w = mh_step_w        # Interval for MH shift proposals
)

# The result 'Res' contains:
# Res$Theta_0   - MAP estimate of the background distribution
# Res$Theta     - List of MAP motif PWMs (probability matrices)
# Res$W         - MAP motif type assignments for each sequence
# Res$A         - MAP motif positions in each sequence
# Res$plot      - A list of ggplot2 objects (log-likelihood trace and logo plots)

# -------------------------------
# Predict motif positions using posterior samples
# -------------------------------

# Use posterior samples to estimate motif type and position for each sequence
pred_frame = pred_ana(
  Data = sample_data,
  theta_0_samp = Res$Theta_0,       # Background distribution
  theta_samp = Res$Theta,           # Motif PWMs
  Motif_length = Motif_length,      # Motif lengths
  burnin = burnin,                  # Burn-in period
  sum_N = sum_N                     # Total number of samples
)

# View predicted motif type and position
print(pred_frame$est_W)   # Predicted motif class (1 to num_motif, or background)
print(pred_frame$est_A)   # Predicted motif location(s) within each sequence

```

## Reference
-   Tang, X., Liu, R. (2025+), “GAMMA: Gap-aware Motif Mining under Incomplete Labeling with Applications to MHC Motifs,” Working Paper.

## Contact
Xinyi Tang: xytang@link.cuhk.edu.hk; Ran Liu: ranliu@bnu.edu.cn
