# Expected GO annotation file
---

The file must be a dataframe with the gene IDs contained in the first row, and GO IDs in the second one.
The order is important.

Gene IDs might be duplicated when the have several GO IDs.

A format like this is expected :

```
geneid	goid
Lalb_Chr00c01g0403641	GO:0005524
Lalb_Chr00c01g0403641	GO:0006464
Lalb_Chr00c01g0403661	GO:0016021
Lalb_Chr00c01g0403671	GO:0016021
Lalb_Chr00c01g0403691	GO:0016021
Lalb_Chr00c01g0403721	GO:0005634
Lalb_Chr00c01g0403721	GO:0003677
Lalb_Chr00c01g0403721	GO:0018108
Lalb_Chr00c01g0403721	GO:0046872
```
