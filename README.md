# phenoGWAS

## Analysis for special variants in interest
5 models
1. lgsbl: Logistic regression for baseline binomial traits
2. coxhm: Cox hazard model until the patient has the binomial outcome. 
3. lnmxi: Linear mixed effect modeling of snp on intercept for the (continous) outcome
4. lnmxs: Linear mixed effect modeling of snp on slope (time-dependent-change) of the outcome
5. lncns: Conditional linear modeling of snp on slope (essentially similar to lnmxi, but the model is expansion of paired t-test, )
    
In total, 52 models were analyzed for each variants.

## Visualized the strength of the associations
92 variants in Meta-5.    
$n = number\ of\ models\ tested$ = 52


1. $0.05 \le P$
2. $\frac{0.05}{n} \le P < 0.05$
3. $\frac{0.05}{92} \le P < \frac{0.05}{n}$     
4. $\frac{0.05}{92 * n} \le P < \frac{0.05}{92}$
5. $ P < \frac{0.05}{92 * n}$

The interpretation would be,
1. no signal    
2. raw P < 0.05
3. FDR 0.05 across traits
4. FDR 0.05 across 92 variants 
5. FDR 0.05 across traits x 92 variants

![Figure1](fig/output.png)
Figure1: 
Variants on x and models on y. 5 colours indicate the levels of significance defined as above. The same traits analyzed by different models are close each other. So a horizontal line indicates somewhat consistent phenotype-genotype associations.    
A couple of interesting points. conditional linear model for slope is the most sensitive analysis. (returns a most of darker cells.)    
APOE1 is associated with insomnia, MMSE, and UPDRS1 seems sensible. 
    
I also conducted GWAS (imputated integer genotypes filtered by maf>0.05 and Rsq > 0.8) but at this moment, models with binomial outcomes (lgsbl/coxhm) were reporting only sub-significant variants. (5e-7). And lncns analyses reported too many variants (2K+). Probably need to filter them with stricter threshold. 