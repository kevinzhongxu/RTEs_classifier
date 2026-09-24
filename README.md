# RTEs_classifier

**RTEs_classifier** is a lightweight, rule-based R script for classifying *reverse-transcribing elements (RTEs)*, represented by *reverse-transcriptase (RT)-containing contig* sequences, based on functional annotations in a GFF3 file.

The classifier integrates RT subtype annotations and associated protein/domain architecture to distinguish major classes of RTEs while retaining ambiguous or weakly supported sequences as separate categories.

RTEs_classifier is intended as a rule-based screening and classification tool rather than a definitive taxonomic assignment method. Classification confidence depends on sequence completeness and annotation quality, and divergent or incompletely annotated elements may remain ambiguous or unclassified.

## Classification

RT-containing contigs can be assigned to:

-   *Retrotransposon-like*
    -   non-LTR
    -   LTR
-   *Retrovirus-related*
-   *Group II intron-like*
-   *Cellular RT-containing*
    -   telomerase-like
-   *Other retroelements*
    -   retron-like
    -   DGR-like
    -   CRISPR-RT-like
-   *Ambiguous RT-containing*
-   *Unclassified RT-containing*

Classification is based on weighted evidence from annotations such as RT subtype, RNase H, integrase, endonuclease, retroviral protease, Gag/capsid, Env, and other characteristic domains.

## Requirements

-   **R**
-   No additional R packages are required.

## Usage

``` bash
Rscript classify_rt_contigs.R input.gff3 [output_prefix]
```

For example:

``` bash
Rscript classify_rt_contigs.R seq_genomes.gff3 seq_genomes
```

The input should be a standard *9-column GFF3 file* containing functional/domain annotations.

## Output

The script generates three tab-delimited files:

``` text
<prefix>_rt_contig_classification.tsv
<prefix>_rt_feature_evidence.tsv
<prefix>_rt_class_summary.tsv
```

These contain:

-   Contig-level classification and confidence
-   Feature-level evidence supporting each assignment
-   Summary counts for the identified RTEs classes

## Classification strategy

The classifier:

1.  Requires prior functional annotation of the sequences, with a standard *.gff3* file as input.
2.  Requires detectable *reverse transcriptase evidence* in the sequence.
3.  Uses combinations of associated domains as supporting evidence.
4.  Reports ambiguous classifications explicitly rather than forcing a sequence into a single class.

## Citation

Baral A, Kothiala A, Zhong KX, Cheung J, Flibotte S, Srebnik S, Dao Duc K, Suttle CA, and Jan E. (2026) Mapping the Functional Landscape and Viral Diversity of 2A Peptides. bioRxiv. DOI: 10.64898/2026.09.24.754024. <https://www.biorxiv.org/content/10.64898/2026.09.24.754024v1>



