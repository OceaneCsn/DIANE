# Regulatory links thresholding

---

Without thresholding, we would obtain a fully cnnected weighted graph from GENIE3, with far too many links to be interpretable.

In order build a meaningfull network, this weighted adjacency matrix betwen regulators and targets has to be sparsified, and we have to determine the regulatory weights that we consider significant.

A first approach is hard thresholding on a desired number of links. If we want a network with 1000 edges, then we choose the 1000 strongest interactions. Usually, gene networks have have a number of edges comparable to their number of nodes, or a little but more. As a first and naive solution, we propose to thresold the network to a number of edges being 1.5*V with V being the number of nodes.
Of course, you can modify the value to use your own threshold.

A more sophisticated method turning regulatory links into pvalues and allowing for significance testing on those links is currently being implemented and will be proposed soon.