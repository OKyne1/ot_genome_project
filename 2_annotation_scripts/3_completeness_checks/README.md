# Completeness Checks
The aim of this step is to identify complete RAGE genes important for classification. This is all tra genes and integrase (we may also included dnaA)

## Processing Steps
1. Conversion of the input genbank to faa file format
2. Blast search of the faa file against a curated database
3. Processing of the blast output to identify genes with >95% alignment (considered complete)
4. Writing this information to the genbank file

At some point we may also want to write traK-1 and traK-2 ireespecitive of it they are complete or not. But this hasn't been handled yet and may not be important.
