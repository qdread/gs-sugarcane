# Genomic selection of sugarcane

## Project contributors

- Quentin Read, USDA ARS
- Md Sariful Islam, USDA ARS


## Notes

- Old code does not seem to work anymore for ADE. The function `sommer::A.mat()` no longer has the argument `shrink`. It did when the code was written, for example see <https://github.com/covaruber/sommer/blob/0333ae177a2304fdeb01e7f668c3db5da57c520d/R/FUN_relationships.R>. Either roll back to an older version of `sommer` or rewrite this code to get the same results with the new version.
- Code could be optimized possibly by vectorizing. Convert each section of the for loop to a function so it can be run on the HPCC. This could help explore a grid of tuning parameters for the ML models.
