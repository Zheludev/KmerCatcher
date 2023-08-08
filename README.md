# KmerCatcher

A simple python script that counts k-mer matches between reference database and query .fasta files, optionally returning all non-matching query sequences. Designed originally for sequencing adapter decontamination by AZF.

=================================================

## dependencies

`KmerCatcher` requires [`tqdm`](https://github.com/tqdm/tqdm) (though not-really, just remove the `tqdm` command and it should all work fine)

=================================================

## overview

Given two input `.fasta` files, a "database" and a "query", `KmerCatcher` will rapidly look for any perfect k-mer matches between the "database" and the "query".

Illumina sequencing adapters can find their way into sequencing data and aren't always perfectly removed. This leads to odd assemblies that are interspersed with bits of adapters. Because these adapters are largely synthetic, we can conclude that if an adapter k-mer is seen in a sequence, it likely is erroneously there and so that sequence should be discarded.

`KmerCatcher` can be used to identify sequences devoid of Illumina adapters by using the provided `illumina_adapters.fa` as the "database" file.

`illumina_adapters.fa` was manually composed first by AZF using [these sequences](https://wikis.utexas.edu/pages/viewpage.action?pageId=28165137) and then I added all the sequences I could find [here](https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/AdapterSequencesIntro.htm)

e.g.:
```
python KmerCatcher.py -k_len 12 -db illumina_adapters.fa -i query.fa -o out.fa -debug 1
```

options:
`-k_len`  :  length (nt) of k-mer used for matching (default = 12)

`-db`  :  input "database" `.fasta` file

`-i`  :  input "query" `.fasta` file

`-o`  :  output zero-match `.fasta` file (if not specified, zero-matches will not be saved)

`-debug`  :  print all match statistics for each sequence as well as some summary notes (default = 1 = on)
