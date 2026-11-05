# Comparative genomics of host association in small-spored *Alternaria*

A pipeline for comparing fungal isolates from different host plants: assembly, gene prediction, functional annotation (CAZymes, effectors, biosynthetic gene clusters), orthology, and host-enrichment statistics.

> [!IMPORTANT]
> **This repo contains no real data or results.** The original project used unpublished isolate data that belongs to someone else, and everything in `data/dummy/` is simulated. The output of the demo says nothing about biology.

> [!NOTE]
> **Reconstruction.** The original scripts from this MSc course project were lost. This code was rebuilt in November 2026, with AI assistance, from the methods section of the project report. It is not the code that produced the report's results.

## Try it (needs R only)

```bash
bash run_demo.sh
```

This simulates tool outputs for 24 made-up isolates, runs steps 08–10, and checks that the planted host-associated features are recovered. Results go to `demo_output/`.

## Pipeline

| Step | Script | Tools |
|---|---|---|
| 1 | `scripts/01_trim_assemble.sh` | Trimmomatic 0.39, SPAdes 3.15.3 |
| 2 | `scripts/02_assembly_qc.sh` | QUAST 5.2.0, BUSCO 5.4.7 (`ascomycota_odb10`) |
| 3 | `scripts/03_predict_genes.sh` | Funannotate 1.8.1 |
| 4 | `scripts/04_functional_annotation.sh` | eggNOG-mapper 2, InterProScan 5.52 |
| 5 | `scripts/05_secretome_effectors.sh` | SignalP 6.0 (slow mode), EffectorP 3.0 |
| 6 | `scripts/06_cazymes_bgcs.sh` | dbCAN3, antiSMASH 6.1.1 |
| 7 | `scripts/07_orthofinder.sh` | OrthoFinder 2.5.4 (DIAMOND + MCL) |
| 8 | `R/08_build_feature_tables.R` | base R |
| 9 | `R/09_host_enrichment.R` | base R |
| 10 | `R/10_plot_counts.R` | base R |

## Run on your own data

1. Install the tools above. SignalP and GeneMark need their own licences.
2. Copy `config/samples_template.tsv` to `config/samples.tsv` and list your isolates, host groups and reads.
3. Edit the paths in `config/config.sh`.
4. Submit steps 1–7 in order from the repo root, for example `sbatch scripts/01_trim_assemble.sh`.
5. Run steps 8–10 using the commands in each script's header.

## Run on public genomes (CAU cluster)

```bash
bash cluster/run_all.sh    # from the repo root, on a login node, inside tmux
```

This creates one conda environment per tool group, downloads 12 public *Alternaria* runs from ENA (`cluster/public_samples.tsv`) and the databases, then submits the full run as one Slurm job. Finished steps are skipped if you resubmit.

This is a test of the pipeline, not a biological result: the two host groups also differ by species and by study. Step 4 is left out, and effector predictions stay empty unless you install SignalP yourself (it needs a licence).

## Method choices

- **CAZymes:** a family is assigned when at least 2 of the 3 dbCAN tools report it; subfamilies are collapsed to family level.
- **Effectors:** EffectorP predictions among SignalP-positive proteins only. Dual predictions are counted under the first-named class.
- **Enrichment:** Fisher's exact test on presence/absence, each host group against the reference host (default `wild_tomato`), with Benjamini-Hochberg FDR across all tests.
- **Candidates:** raw p < 0.01 and odds ratio > 5. These are exploratory, not significant discoveries.

## What has and hasn't been tested

| Part | Status |
|---|---|
| Steps 8–10 | Run on simulated data, with checks in `tests/check_demo.R` |
| Steps 1–7 | Script logic tested with stand-in programs only; **never run with the real tools** |
| `cluster/run_all.sh` | Written for the CAU cluster; not yet run there |
| Parsers in step 8 | Written from the documented output formats and not yet tried on real tool output, so check them against your tool versions |

## Differences from the report

- SPAdes has no `--fungal` option, so step 1 uses `--careful` alone.
- The report doesn't give these, so they are assumptions: the Trimmomatic `ILLUMINACLIP` settings and the Funannotate `clean`, `sort` and `mask` steps.
