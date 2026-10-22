# Publish this adapter

This standalone variation was derived from the local FSHR example at commit
`1867928`. It uses *HFE*, a different gene involved in iron metabolism, and its own
Ensembl data rather than relabeling the FSHR data.

Gene metadata and variant responses were downloaded on 2026-10-22 from:

- [HFE gene lookup](https://rest.ensembl.org/lookup/symbol/homo_sapiens/HFE?content-type=application/json)
- [rs1800562 variation data](https://rest.ensembl.org/variation/human/rs1800562?pops=1&genotypes=1&content-type=application/json)
- [rs1799945 variation data](https://rest.ensembl.org/variation/human/rs1799945?pops=1&genotypes=1&content-type=application/json)

The response snapshots are in `data/`. Gene coordinates use GRCh38. The graph
schema, adapter interface, tests, and Croissant record sets follow the FSHR
example's format. Cached data and graph outputs are included by `.gitignore`.

New and adapted functions carry `AI-Generated` comments for manual review.

TODO: replace placeholder author and GitHub URLs in `pyproject.toml` with the
actual publication details.

TODO: confirm the source code license and upstream data terms. The original
example's README said MIT, but no LICENSE file was supplied. No dataset license,
publisher, publication date, or repository URL has been invented in Croissant.

TODO: add the preferred citation, published repository URL, confirmed data
license, and publication date to `export_croissant.py` when those are available.
