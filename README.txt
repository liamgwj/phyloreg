
NEEDS UPDATING

Current workflow:

home directory: 'phyloreg'

folder 'datasets' contains subfolders for various data sources
each subfolder contains:
1. a host/pest compatibility matrix (referred to as 'database' or 'db')
2. a host phylogeny
3. cleaning scripts and associated files used to produce items 1 & 2

data format:
'db' is a binary matrix where rows represent pest/symbiont species and columns
represent host species. '1's indicate that a pair are capable of associating
while '0's indicate (assumed) non-compatibility. Rows and columns are named,
and all rows and all columns sum to at least 1 (i.e. every pest has at least
one known host and vice versa).
'phy' is the host phylogeny, where all and only host species from 'db' are
present as a tips.


folder 'scripts' contains generic scripts that are run from source by the '_UI'
r scripts in the top-level directory


'model_UI.r' runs R script 'scripts/logreg_generic.r', which fits logistic
regression models for all host/pest combinations in 'db' and writes output to
'output/model_coefficients'.

'plotting_UI.r' runs 


Premise: ----------------------------------------------------------------------

Observation: using logistic regression to predict the probability of two hosts
sharing a pathogen given the PD separating them relies on '0' observations
being true non-associations

In the case where '0's are assumed rather than measured, the model predictions
are likely to be influenced by the number and/or distribution of '0'
"observations" that are included

We want to assess the potential variation in model predictions under different
false-zero situations



Method:

start with Robles (code available) then Gilbert (bigger reach)


1. reproduce published analysis (logistic regression portion only). This will
include some number of '0' observations.

2. re-run analysis with successively fewer, then more '0' observations, writing
out slopes/predictions at each step.

3. compare outputs of varied models

