# Gene information file requirements

---
This file gives additional information about the genes in your analysis.

It must contain a column named "Gene".
This column should contain the genes IDs present in the row names of the expression file, or the locus corresponding to your splice variants. Genes can be duplicated, in that case the final annotation fields will be a concatenation of the common annotations.

The annotation fields must be either "label", "description", or both.
The column "label" is typically used for common gene names. "description" is usually a brief description of the known roles/pathways of the gene.

The summary tables containing genes will contain those additional columns.

Here is an example of a tab separated file that would work :
```
label   description Gene
HTT2	target of trans acting-siR480/255 protein	AT5G18040
HRS1	histidyl-tRNA synthetase	AT3G46100
HIPP20	Heavy metal transport/detoxification superfamily protein	AT1G71050
```