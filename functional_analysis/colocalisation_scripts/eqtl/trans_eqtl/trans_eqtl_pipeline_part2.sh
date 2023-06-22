# Has to be run in the base directory for the analysis, and results after
# sending away part 1 output need to be in part2_input/
#
# bash trans_eqtl_pipeline_part2.sh

Rscript tqtl_coloc.R $(pwd)

bash annotate_ensembl_gene.sh
echo "Part 2 done - check Coloc_results/ for results"

