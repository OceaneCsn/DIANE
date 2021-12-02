# Version changes in DIANE
---

### 1.0.1 : DIANE is released

Version at the moment of the official BMC publication (May 26, 2021). 


### 1.0.2 : Minor bug fixes and improvements.

Improvements :

+ Citation was added
+ Gene IDs can now have blank space round them in the exploration of specific genes profiles
+ MDS was removed
+ PCA can now handle expression data with 4 or less experimental conditions
+ This document was added to the interface in a new tab

Bug fixes : 

+ Custom annotation/gene information file now allows duplicated genes
+ Custom annotation/gene information now can only have columns named label and/or description
+ DEA csv download does not require annotations any more

### 1.0.3 : Minor bug fixes and improvements.

+ Fixed error in DAINE's programming interface (wrong value of the defaut "conditions"" argument in normalize function) 22/11/2021
+ Used namesapce operator in data import module to speed up loading