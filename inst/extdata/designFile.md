# Expected design file

---

You can also feed DIANE the experimental design corresponding to your data.

Indeed, DIANE can anlysis any combinatorial design and give valuable insights about the different factors that are studied.

To do so, you can provide, for each condition name, the level of the factors corresponding to that condition in your study.


Let's take our companion data as an example.
Here, the experiment we chose for our demo included 3 factors, heat stress, mannitol stress, as well as salinity stress.

 
You should first identify which level of each factor can be considered as the **control** level, and which is a **perturbation**. In our demo, the control condition is referred to as C.
Then, the perturbation are specified by their corresponding letter, from simple to triple stress combination. For example, SH corresponds to salt and heat stresses in the control level of mannitol. As a consequence, its levels would be 1,0,1.

The design file is thus a matrix with condition names as rows, and factor names as columns.
It MUST contain a column named "Condition", as you can see in our example below :

```
Condition,Salt_stress,Mannitol_stress,Heat_stress
C,0,0,0
H,0,0,1
S,1,0,0
M,0,1,0
SM,1,1,0
SH,1,0,1
MH,0,1,1
SMH,1,1,1
```
Many RNA-Seq anlysis only study one factor, which would not be a problem to use DIANE. If you have two conditions named control and trt, the design would contain only one factor column: 

```
Condition,treatment
contr,0
trt,1
```