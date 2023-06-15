# Altered script from mkanai/slalom-paper repo

import hail as hl

for files in ["data/cup_file/FASTA_BED.ALL_GRCh37.novel_CUPs.bed",
              "data/cup_file/FASTA_BED.ALL_GRCh37.reject_2.bed"]:
    ht = files.replace(".bed", ".ht")
    bed = hl.import_bed(files, reference_genome="GRCh37", min_partitions=10)
    bed.write(ht)

for files in ["data/cup_file/FASTA_BED.ALL_GRCh38.novel_CUPs.bed",
              "data/cup_file/FASTA_BED.ALL_GRCh38.reject_2.bed"]:
    ht = files.replace(".bed", ".ht")
    bed = hl.import_bed(files, reference_genome="GRCh38", min_partitions=10)
    bed.write(ht)
