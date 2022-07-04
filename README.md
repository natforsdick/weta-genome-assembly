# Wētā genome assembly

All scripts related to genome assembly for Little Barrier giant wētā (*Deinacrida heteracantha*).

Input data: PacBio CLR 'HiFi-like', Illumina HiSeq, Hi-C.

[01-assembly](01-assembly/) - testing hifiasm, HiCanu, and MaSuRCA to assemble PacBio 'HiFi-like' data. To date (2022-07-04) only hifiasm has completed.

[02-purge-dups](02-purge-dups/) - Implementing the purge-dups pipeline to remove duplicate contigs.

[03-polish](03-polish/) - Using HiSeq data to polish the purged assembly.

[04-hic-arima](04-hic-arima/) - Using the Arima pipeline to map Hi-C data to the purged assembly, then SALSA or YAHS to perform scaffolding.

[05-alignment](05-alignment/) - Aligning draft genomes against one another for comparisons.

[06-merqury](06-merqury/) - Implementing Merqury pipeline for assembly QC.

[07-read-correction](07-read-correction) - Attempting correction of raw PacBio CLR reads using FMLRC2, with the intention of using a draft assembly as a 'scaffold' on which to assemble these reads against. 
