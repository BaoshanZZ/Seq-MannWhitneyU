# R Scripts for Publication:  
**"Group Sequential Test for Two-Sample Ordinal Outcome Measures"**

This folder includes the following major R scripts for our publication.

## R Scripts Overview

### **Func_Tests.R**  
This R script includes functions for different test statistics and for the t-adjusted critical boundary.

#### **Functions:**
- `calc.Zknew()`: For our proposed test  
- `wilcoxon.h0()`: For Wilcoxon test with variance evaluated assuming two arms have the same distribution  
- `brunner()`: For Brunner and Munzel test  
- `brunner.log.ratio()`: For Brunner and Munzel test with log odds transformation  
- `ctz.adj()`: For t-adjusted critical boundary  

---

### **Sim_Category.R**  
This R script conducts the simulation for categorical data, comparing different test statistics by calling corresponding functions from `Func_Tests.R`.

---

### **Sim_Exponential.R**  
This R script conducts the simulation for skewed continuous data (with an exponential distribution), comparing different test statistics by calling corresponding functions from `Func_Tests.R`.

---

### **Sim_Poisson.R**  
This R script conducts the simulation for count data (with a Poisson distribution), comparing different test statistics by calling corresponding functions from `Func_Tests.R`.
