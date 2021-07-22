# Codes for analyzing the SARS-CoV-2 Delta variant outbreak on May 21 2021 in Guangzhou, Guangdong, China
1. The sequencing analysis pipeline (variant_pipeline.py) was used for calling variants and generated consensus sequences from the raw sequencing data from illumina; We used the tools in https://github.com/ItokawaK/Alt_nCov2019_primers to trim primers and iVar https://github.com/andersen-lab/ivar to call variants and generate the consensus sequence.
2. Tree was bulid with phyml, and annotated by using R package ggtree (R _script_tree.R).
