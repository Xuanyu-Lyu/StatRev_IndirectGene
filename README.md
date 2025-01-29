# StatRev_IndirectGene
A repository for the statistical review of all the indirect genetic effects models. 

## Research Goal
Simulate indirect genetic effects with a modified version of GeneEvolve. Use the simulated data to fit as many models on indirect genetic effects as possible to compare the bias and accuracy of the estimates under various conditions (including but not limited to sample sizes, the strength of PGS and the effect sizes of different indirect genetic effects). Come up with a clear guideline for future empirical research on indirect genetic effects from a modeling standpoint (in contrast to a phenotypic standpoint).

## Target Journal
I would like to submit the paper to a higher-impact journal, but it is not too far for me. For example, **Nature Communications, Science Advance, Nature Human Behavior, Annual Review of Genomics and Human Genetics** will be good places for this paper. I'm open to other suggestions, so please contact me if there are other better options in mind.

## Involved Collaborators
Not ranked.
- Xuanyu Lyu: Design Study, Modify Simulation Functions, Run the Simulations, Fit models, Write Manuscripts
- Tong: Modify Simulation Functions, Help with model fitting and interpretations
- Dinka: Design Study, Propose Models, Interpret Results
- Noemie: Design Study, Propose Models, Interpret Results
- Richard Border: Help with simulations (depend on what types of effects we want to simulate/ also we need to go back to the bias from his simulations)
- Matt: Captain Responsibilities

## Project Specifics
### What models should we include?
- PGS-Based Models
  - SEM-PGS
  - Regression-based Family PGS approaches:
    - Y ~ PGS_o + PGS_p + ...
    - Y ~ PGS_within + PGS_between + ...
    - Y ~ PGS_NT + ...
    - Y ~ PGS_adoptee + PGS_Biological + ...
  - MR-DoC2
- Classical Family Designs
  - Children of Twins
  - Cascade
- Models that want to control for IGE (Indirect genetic effects)
  - Trio-GCTA
  - RDR
  - With-Sib GWAS
- Prenatal Indirect Genetic effects (need more wisdom here!)

### What types of indirect genetic effects do we want to simulate?
Bottom line: we want to simulate multivariate effects, no matter what kind of effects we are going to simulate
- Genetic Nurture as a consequence of vertical transmission- different vertical transmission from mother and father. e.g., Prenatal Effects
- Sibling Nurture
- Horizontal Pleitropy (?) Go back to the issue that David Evans
- Assortative Mating
- Environmental Effects Span Generations
- GE interaction (possible)
- G-E covariance caused by the active choice of parents that manifests in both parent and offspring generation

### What coding language should we use to simulate?
Python or R. One other goal I'm anticipating for the project is to publish a package based on Gene-Evolve that is dedicated to simulating different types of indirect genetic effects. 
- Pros for Python: Faster, Personal Growth for Xuanyu, Easier to Migrate to other projects (specifically the prediction neural network project I'm planning on).
- Cons for Python: Limited Users in BG. We still probably need OpenMx to fit most of the model, so it would be good if the model fitting could be included in the same package in R.
- Pros for R: We have a lot of codes available. It is easier for me and tong to implement any changes.
- Cons for R: Slower.
Matt's Opinion: Use Richard's package if he's available to implement simulations on other indirect genetic effects that have not been covered yet. 
