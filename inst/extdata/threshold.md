# Regulatory links thresholding

---

Without thresholding, we would obtain a fully cnnected weighted graph from GENIE3, with far too many links to be interpretable.

In order build a meaningfull network, this weighted adjacency matrix betwen regulators and targets has to be sparsified, and we have to determine the regulatory weights that we consider significant.