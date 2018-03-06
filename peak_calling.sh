python ../../coverage_mean.py ../$control"_peaks.bed" $control"_av_50nt.bed"
sort -k1,1 -k2,2n $control"_av_50nt.bed" > $control"_av_srt_50nt.bed"
rm $control"_av_50nt.bed"
bedtools intersect -wo -a ../$treatment"_peaks.bed" -b $control"_av_srt_50nt.bed" > $experiment"_1vs50_wo_is.bed"
cat $experiment"_1vs50_wo_is.bed"\
 | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,($4-$8>0)?$4-$8:-($4-$8)}'\
 >$experiment"_overlap_1vs50_wo_is.bed"
 cat $experiment"_overlap_1vs50_wo_is.bed"\
  | awk 'BEGIN{FS=OFS="\t"} {if($9>="'${thr}'"+0){print $1,$2,$3,$9}}'\
  > $experiment"_thr"$thr"_peaks.bed"
 bedtools intersect -v -a ../$treatment"_peaks.bed" -b $control"_av_srt_50nt.bed" > $experiment"_nonov_peaks.bed"
 cat $experiment"_nonov_peaks.bed" | awk '{if($4>="'${thr}'"+0){print}}'\
  >> $experiment"_thr"$thr"_peaks.bed"
 sort -k1,1 -k2,2n $experiment"_thr"$thr"_peaks.bed" > $experiment"_srt_thr"$thr"_peaks.bed"
