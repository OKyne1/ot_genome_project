# MLST
MLST (Multilocus Sequence Typing) is a molecular typing method used to characterize bacterial isolates by analyzing the sequences of several (usually 6-7) housekeeping genes. This method provides a way to classify strains based on their unique sequence types. To date all _Orientia_ genomes have had different combinations of MLST genes.

For assemblies MLST can be used to check for genome completeness by confirming that all targeted housekeeping genes are present and intact in the sequenced genome. This is complementary to BUSCO as they use different genes.

Usage:
`mlst ../genomes/*.fasta > mlst_output.txt`

The output will look something like this:
wgot001.fasta        otsutsugamushi  -       gpsA(26)        mdh(1)  nrdB(8) nuoF(39)        ppdK(29)        sucB(22)             sucD(~35)

wgot016.fasta        otsutsugamushi  -       gpsA(26)        mdh(1)  nrdB(3) nuoF(10)        ppdK(7) sucB(24)        sucD(22)

wgot031.fasta        otsutsugamushi  -       gpsA(36)        mdh(1)  nrdB(3) nuoF(6) ppdK(36)        sucB(24)        sucD(22)

Tilde symbols mean that there is potentially a new gene allele. New genomes should probably have their genes submitted to MLST, talk to Liz about this.