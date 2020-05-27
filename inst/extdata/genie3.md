
# GENIE3 network inference

---

The [GENIE3 package](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0012776) is a method mased on machine learning to infer regulatory links between genes and regulators.

GENIE3 needs to be fed a list of genes, that will be the nodes of the inferred network. Among those genes, some must be considered as potential regulators. 

GENIE3 can determine the influence if every regulators over each input genes, using their respective expression profiles. You can specify which conditions you want to be consired for those profiles during the network inference.

For each target gene, the methods uses Random Forests to provide a ranking of all regulators based on their influence on the target expression. This ranking is then merged across all targets, giving a global regulatory links ranking.

The idea is then to keep the strongest links to build the gene regulatory network. The way of choosing this minimal importance value needed to be included in the network will be described the "thresholding" box.

GENIE3 was among the best performers of the [DREAM challenge](https://www.synapse.org/#!Synapse:syn2787209/wiki/70350), designed to benchmarck state of the art network inference methods of simulated and validated biological data. The advantages of the method is that it only gives oriented edges from regulators to targets, which is desired in the context of regultory networks, and captures regulators interactions and combinatorics. 