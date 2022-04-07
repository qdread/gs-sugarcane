# Genomic selection of sugarcane

## Project contributors

- Quentin Read, USDA ARS/NC State University
- Md Sariful Islam, USDA ARS
- Alex Lipka, University of Illinois
- Per McCord, Washington State University

## How to run

This requires phenotype and genotype data provided by Sarif. The main script is `GS_sugarcane_2022.R` which needs to be run on a remote cluster. Functions called in that script are sourced from `GS_sugarcane_2022_fns.R`. Additional scripts and notebooks are for exporting and visualizing output.

## Notes

- Old code does not seem to work anymore for ADE. The function `sommer::A.mat()` no longer has the argument `shrink`. It did when the code was written, for example see <https://github.com/covaruber/sommer/blob/0333ae177a2304fdeb01e7f668c3db5da57c520d/R/FUN_relationships.R>. Either roll back to an older version of `sommer` or rewrite this code to get the same results with the new version.

