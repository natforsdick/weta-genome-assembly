# Wētā genome assembly

All scripts related to genome assembly for Little Barrier giant wētā (*Deinacrida heteracantha*). This project is a Genomics Aotearoa & MWLR collaboration led by Thomas Buckley. Contributors include Manpreet Dhami, Ann Mc Cartney, Dukchul Park, Dini Senanayake, Natalie Forsdick.

This repo doesn't contain anything particularly novel, it merely is a record of the various tools used for wētā genome assembly. No doubt there are improvements that could be made to my scripting.

Input data: PacBio CLR 'HiFi-like', Illumina HiSeq, Hi-C.

[01-assembly](01-assembly/) - Testing hifiasm, HiCanu, and MaSuRCA to assemble PacBio 'HiFi-like' data. To date (2022-07-04) only hifiasm has completed.

[02-purge-dups](02-purge-dups/) - Implementing the [purge-dups pipeline](https://github.com/dfguan/purge_dups) to remove duplicate contigs.

[03-scaffolding](03-scaffolding/) - Using the [Dovetail Omni-C pipeline](https://omni-c.readthedocs.io/en/latest/index.html) and [YaHS](https://github.com/c-zhou/yahs) to scaffold the assembly.

Subsequent steps are all just exploratory at this point.

[04-fill-polish](04-fill-polish/) - Exploring the use of gapfilling tools to improve the scaffolded assembly.

[05-alignment](05-alignment/) - Aligning draft genomes against one another for comparisons.

[06-read-correction](06-read-correction) - Attempting correction of raw PacBio CLR reads using FMLRC2, with the intention of using a draft assembly as a 'scaffold' on which to assemble these reads against. 

[QC](QC/) - Contains scripts for raw read QC and assembly QC - fastqc, Hi-C and Pore-C QC, BUSCO5, Merqury.
