# Variable grouping in network inference

---

As the inference algorithm determines the influence of each regulator on each target, high correlation between regulators, or variables, can be problematic.
For instance, only one of the two very correlated regulators could be detected by the random forests, losing the information of the other interaction.
Or, their influence could be shared, resulting in a spuriously low importance value for the two pairs.

To prevent that, we consired grouping the regulators with a Spearman (non linear correlation) exceeding a certain value, 90% as default.

New regulators are created, called grouped regulators, and will be represented in the graph as darker and bigger square nodes.