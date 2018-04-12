## SICER

When using SICER without control, you should run SICER-rb.sh which requires the strand information. As BLISS data has no strand info, you should put some to it.
`cat $ID"_exp.bed" | awk '{ print $0 "\t0\t+"}' > $ID"_exp_full.bed"`

Than you run SICER-rb.sh
`$sicerpath"SICER-rb.sh" $inpath BB101_Milan-1-77_exp_full.bed $outpath hg19 100 200 150 0.8 600 0.1 >$log`

Converting wig to bed with [bedops](https://bedops.readthedocs.io/en/latest/)
`wig2bed < $ID"_exp_full-W200-G600-E0.1-islandfiltered-normalized.wig" > $ID"exp_full-W200-G600-E0.1-islandfiltered-normalized.bed"`
