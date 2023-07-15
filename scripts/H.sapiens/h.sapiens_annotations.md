# Annotations for Homo Sapiens

## rRNA
rRNA annotations were obtained from ITAS

Stupnikov, A., Bezuglov, V., Skakov, I., Shtratnikova, V., Pilsner, J. R., Suvorov, A., & Sergeyev, O. (2022). ITAS: Integrated Transcript Annotation for Small RNA. Non-coding RNA, 8(3), 30. https://doi.org/10.3390/ncrna8030030

Filtered for rRNA.

## snoRNA
snoRNA annotations were obtained from snoDB.

Bouchard-Bourelle P, Desjardins-Henri C, Mathurin-St-Pierre D, Deschamps-Francoeur G, Fafard-Couture É, Garant JM, et al. snoDB: an interactive database of human snoRNA sequences, abundance and interactions. Nucleic Acids Research. 2020 Jan 8;48(D1):D220–5. 

Converted using, 

```
python conversion_tools.py tsv_to_gtf ...
```

This command prompts the user to fill in information which makes it easy to recreate. 

One can even get box specific information by combining box info

```
python gtf_modifiers.py merge_two_attributes ...
```

Please consult documentation for proper use. 

## snRNA
snRNA annotations were obtained from RNA central.

Filtered for transcripts and snRNAs. 

RNAcentral Consortium, Sweeney BA, Petrov AI, Ribas CE, Finn RD, Bateman A, et al. RNAcentral 2021: secondary structure integration, improved sequence search and new member databases. Nucleic Acids Research. 2021 Jan 8;49(D1):D212–20. 
