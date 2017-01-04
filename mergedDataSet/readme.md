# BRC_SkinGut_Microbiome/mergedDataSet
Analysis of the merged data set for a set of questions.

**Project.id = 4**


Scripts in no particular order.  

1. **import.R**  
..* Import the 3 month old data with covariates, look at the distributions of diversity. Look at the relationships of diversity density distributions with eczema status at 3 months conditional on covariates.  
2. **univariate_report.R**  
..* Import the 3 month old data with covariates. Regress the Shannon diversity index and logit transformed abundances of microbes on each of the covariates and report p-values. Split the Eczema 3 month group into low and high diversity and plot the distributions.  
3. **eczema12m_diversity.R**  
..* Model the eczema diversity at 12 months, and look at the associations of any covariates in a logistic regression model. Model the change from 3 months eczema status to 12 months as a multinomial distribution and binomial distribution and report the significant covariates.  
4. **diversity_high_group.R**  
..* Testing for covariates that show significant relationship with diversity high and low groups - nothing really important here.  
5. **diversity_report.R**  
..* multiple regression analysis for diversity at 3 months.  

