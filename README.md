# hfe-variant-prevalence

A BioCypher adapter for population prevalence of the *HFE* variants rs1800562
(p.Cys282Tyr) and rs1799945 (p.His63Asp), using 1000 Genomes data from Ensembl.
*HFE* encodes the homeostatic iron regulator. This is a standalone example using
the same graph model and Croissant format as the FSHR adapter.

## Included data

- `data/ensembl_rs1800562.json` and `data/ensembl_rs1799945.json`: real cached
  Ensembl variation responses, including population frequencies and genotypes.
- `data/ensembl_HFE_gene.json`: the Ensembl gene lookup used for gene metadata.
- `data/genotype_frequencies.csv`: genotype prevalence for each variant.
- `data/variant_combinations.csv`: estimated joint diplotype prevalence.
- `biocypher-out/`: Neo4j-compatible node and edge CSVs with separate headers.
- `croissant.json`: Croissant 1.0 metadata and checksums for the bundled data.

Aggregate ALL genotype frequencies use observed sample counts. Subpopulation
genotype frequencies use Hardy-Weinberg equilibrium estimates from allele
frequencies. Joint diplotype frequencies multiply the per-variant frequencies
under an independence assumption; they are estimates rather than observed joint
frequencies. Use `estimation_method` to distinguish them.

The bundled graph contains 147 genotype-frequency records and 162 combination
records. The original parser estimates subpopulation genotypes only when two
allele-frequency entries are available; monomorphic and non-biallelic population
records are skipped. Consequently the exported graph covers 31 populations and
aggregates, while the raw snapshots retain all 32 phase 3 population groups.

## Run

```bash
uv sync --group dev --extra croissant
uv run pytest tests/ -q
uv run python export_croissant.py
uv run python create_knowledge_graph.py
```

`export_croissant.py` rebuilds both standard CSV tables and their metadata using
only cached responses. The graph builder uses the cached data too, but BioCypher
may download the configured Biolink ontology. Run from this project directory.

The adapter yields `(node_id, label, properties)` nodes and
`(edge_id, source_id, target_id, label, properties)` edges. It supports injected
`raw_data` for deterministic offline tests.

## Graph model

The node types are `gene`, `sequence variant`, `population of individual
organisms`, `genotype frequency`, and `variant combination frequency`.

The edges connect genes to variants, variants to genotype frequencies, genotype
frequencies to populations, combinations to populations, and combinations to
their variants. Both frequency node types expose `prevalence_percent`.

## Croissant

The metadata follows the [Croissant 1.0 specification](https://docs.mlcommons.org/croissant/docs/croissant-spec.html).
Load it locally with the optional Croissant dependency installed:

```python
import mlcroissant as mlc

dataset = mlc.Dataset("croissant.json")
rows = dataset.records(record_set="genotype_frequencies")
```

The record sets use standard comma-delimited CSV. The Neo4j import files use the
same semicolon delimiter and single-quote character as the original FSHR example.

See [PUBLISHING.md](PUBLISHING.md) for provenance and publication details.
