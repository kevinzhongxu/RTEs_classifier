# RTEs_classifier
RTEs_classifier is a lightweight, rule-based R script for classifying reverse-transcriptase (RT)-containing contigs based on functional annotations in a GFF3 file.

The classifier integrates RT subtype annotations and associated protein/domain architecture to distinguish major classes of retroelements while retaining ambiguous or weakly supported sequences as separate categories.

Classification

RT-containing contigs can be assigned to:

Retrotransposon-like
non-LTR
LTR
Retrovirus-related
Group II intron-like
Cellular RTase-containing
telomerase-like
Other retroelements
retron-like
DGR-like
CRISPR-RT-like
Ambiguous RT-containing
Unclassified RT-containing

Classification is based on weighted evidence from annotations such as RT subtype, RNase H, integrase, endonuclease, retroviral protease, Gag/capsid, Env, and other characteristic domains.

Requirements
R
No additional R packages are required.

Usage
Rscript classify_rt_contigs.R input.gff3 [output_prefix]

For example:

Rscript classify_rt_contigs.R seq_genomes.gff3 seq_genomes

The input should be a standard 9-column GFF3 file containing functional/domain annotations.

Output

The script generates three tab-delimited files:

<prefix>_rt_contig_classification.tsv
<prefix>_rt_feature_evidence.tsv
<prefix>_rt_class_summary.tsv

These contain:

Contig-level classification and confidence
Feature-level evidence supporting each assignment
Summary counts for the identified retroelement classes
Classification strategy

The classifier:

Requires detectable reverse transcriptase evidence.
Gives higher weight to element-specific annotations.
Uses combinations of associated domains as supporting evidence.
Reports ambiguous classifications explicitly rather than forcing a sequence into a single class.


Citation


