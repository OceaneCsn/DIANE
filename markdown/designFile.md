# Expected design file

---

You should also feed DIANE the experimental design corresponding to your data.

Indeed, DIANE can anlysis any combinatorial design and give valuable insights about the different factors that are studied.

To do so, you must provide, for each condition name (cNf for example), the level of the factors of you study corresponding to that condition.


Let's take our demo data as an example.
 Here, lower or upper case in condition names stand for the amount of environment factor that the plant was exposed to. We studied 3 factors, atmospheric CO2, Nitrate concentration in the culture medium, and Iron starvation.
 
 For us, "c" means ambiant CO2, and "C" is elevated CO2, same for "n" and "N" representing low or high Nitrate concentrations, and "f" and "F", meaning respectively Iron starvation and iron supply.

 
You should first identify which level of each factor can be considered as the **control** level, and which is a **perturbation**. In our study, ambiant is the control level for CO2 exposure, whereas elevated is the perturbation level. For nutritionnal levels, low amounts of nitrate and iron are perturbations.

Control levels have a value of 0, and perturbations have a value of 1.
Hence, if you take the condition "cNf", the levels of our 3 factors are 0 for CO2, 0 for Nitrate and 1 for the Iron.

The design file is thus a matrix with condition names as rows, and factor names as columns.
It should contain a column named "Condition", as you can see in our example below :

```
Condition,CO2,Nitrate,Iron
cnf,0,1,1
cnF,0,1,0
cNF,0,0,0
cNf,0,0,1
Cnf,1,1,1
CNF,1,0,0
CNf,1,0,1
CnF,1,1,0
```
Many RNASeq anlysis only study one factor, which would not be a problem tu use DIANE. If you have two conditions named contr and trt, the design would contain only one factor column: 

```
Condition,treatment
contr,0
trt,1
```

Please, for now, only comma seperated files are accepted.
