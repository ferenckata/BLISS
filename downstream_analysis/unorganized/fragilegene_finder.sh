 # the output of venn maker misses the enters at the end of each title line, that should be removed
 for file in BB11*fragile*gene.bed;\
 do name=$(echo $file|cut -d"_" -f1);echo $name;\
 cat $file | tr ";" " " |awk '{print $0 "\t" $3-$2}'>$name"_fragile_gene_withlength.tsv";\
 done
