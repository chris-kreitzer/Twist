#!/bin/bash


# misc functions:

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# converting .bcf to .vcf format
bcftools view <input.bcf> | bcftools view > <output.vcf>  # vcf requires much more space than binary bcf

# grep'in NV2.10864 (twist) from .vcf file
grep 'NV2.10864' <search.file.vcf>   

# making symbolic link to data repository
ln -s /proj/ferrer/rna[...] <destinationfolder>

