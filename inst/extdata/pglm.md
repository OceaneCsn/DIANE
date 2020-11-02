# Generalized Poisson regression on a cluster of genes

------------------------------------------------------------------------

The idea is to extract the importance and effect of each factor. To do so, the expression of each gene is modeled by a Poisson distribution. The log of its parameter (the expected value) is approximated by a linear combination of the factors in the experiment.

The coefficients associated to each factors are estimated to fit gene expression, and can be insightful to characterize genes behavior in a particular cluster. The model with interactions is considered automatically. If your design in not a complete crossed design, the interaction term will be null.

This should be interpreted only when a cluster is homogenous.
