
cut -f 1,2,3,4,5 elephant_exon_positions.txt > elephant_all_exons.bed

grep 'protein_coding' elephant_exon_positions.txt > elephant_coding_exons.bed
grep -v 'protein_coding' elephant_exon_positions.txt > elephant_nocoding_exons.bed
cut -f 1,2,3,4,5 elephant_coding_exons.bed > elephant_coding_exons.bed2
cut -f 1,2,3,4,5 elephant_nocoding_exons.bed > elephant_nocoding_exons.bed2

mv elephant_coding_exons.bed2 elephant_coding_exons.bed
mv elephant_nocoding_exons.bed2 elephant_nocoding_exons.bed


cat elephant_restriction_fragments.bed elephant_restriction_fragments.bed2 > elephant_all_restriction_fragments.bed
sort -u elephant_all_restriction_fragments.bed > elephant_all_restriction_fragments.bed2
mv elephant_all_restriction_fragments.bed2 elephant_all_restriction_fragments.bed
sed 's/:/\t/g' elephant_all_restriction_fragments.bed > ID
cut -f 1,2,3 ID > ID2
paste ID2 elephant_all_restriction_fragments.bed > test
mv test elephant_all_restriction_fragments.bed
rm ID ID2
rm elephant_restriction_fragments.bed elephant_restriction_fragments.bed2

cut -f 1 elephant_all_restriction_fragments.bed > chr
cut -f 2,3,4 elephant_all_restriction_fragments.bed > other
sed -i 's/chr//g' chr
paste chr other > elephant_all_restriction_fragments.bed
rm chr other

./overlap_converted_frag.py elephant_all_restriction_fragments.bed elephant_all_exons.bed elephant_overlap_all_exons.txt T

ID.human\tID.opossum\tNbExcludedBases1\tNbExcludedBases2\tTotalUngappedLength\tTotalIdenticalLength\tFilteredUngappedLength\tFilteredIdenticalLength\tTotalAlignmentLength\tFilteredAlignmentLength


cat AlignmentStatistics_TBA_ExcludingGenes.txt_x* > AlignmentStatistics_PECAN_ExcludingGenes.txt
grep -v 'TotalUngappedLength' AlignmentStatistics_PECAN_ExcludingGenes.txt > AlignmentStatistics_PECAN_ExcludingGenes.txt2
sed -i '1i ID.mouse\tID.human\tNbExcludedBases1\tNbExcludedBases2\tTotalUngappedLength\tTotalIdenticalLength\tFilteredUngappedLength\tFilteredIdenticalLength\tTotalAlignmentLength\tFilteredAlignmentLength' AlignmentStatistics_PECAN_ExcludingGenes.txt2
mv AlignmentStatistics_PECAN_ExcludingGenes.txt2 AlignmentStatistics_PECAN_ExcludingGenes.txt
rm AlignmentStatistics_TBA_ExcludingGenes.txt_x*
  