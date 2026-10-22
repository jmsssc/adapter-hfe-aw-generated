#!/usr/bin/env python3
"""
hfe-variant-prevalence - BioCypher knowledge graph builder.

Builds a knowledge graph of the HFE gene, its variants rs1800562 and rs1799945,
their genotype (zygosity) prevalence in 1000 Genomes populations, and the
estimated joint diplotype ("combination") prevalence across both variants.

Data is retrieved from the Ensembl REST API (which serves 1000 Genomes
frequencies) and cached under ``data/``.
"""

import logging

from biocypher import BioCypher

from hfe_variant_prevalence.adapters.hfe_variant_adapter import HfeVariantAdapter

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def main() -> None:
    """
    AI-Generated.
    Main function to create the knowledge graph.
    """
    logger.info("Starting hfe-variant-prevalence knowledge graph creation")

    bc = BioCypher(
        biocypher_config_path="config/biocypher_config.yaml",
        schema_config_path="config/schema_config.yaml",
    )

    # Fetches rs1800562 and rs1799945 from Ensembl (1000 Genomes) and caches to data/.
    adapter = HfeVariantAdapter(
        variants=["rs1800562", "rs1799945"],
        cache_dir="data",
        use_cache=True,
        estimate_combinations=True,
    )

    logger.info("Creating knowledge graph...")
    bc.write_nodes(adapter.get_nodes())
    bc.write_edges(adapter.get_edges())

    logger.info("Knowledge graph creation completed successfully!")
    bc.summary()


if __name__ == "__main__":
    main()
