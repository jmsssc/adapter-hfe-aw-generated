"""
HfeVariantAdapter

BioCypher adapter that models the HFE gene, its variants (default: rs1800562 and
rs1799945), and their population prevalence in the 1000 Genomes Project, retrieved
from the Ensembl REST API.

The adapter produces the following nodes:

* ``gene``                          -> HFE
* ``sequence variant``              -> rs1800562, rs1799945 (Biolink: SequenceVariant)
* ``population of individual organisms`` -> 1000 Genomes populations
* ``genotype frequency`` (leaf)     -> per variant x population x genotype
                                       (homozygous ref / heterozygous /
                                       homozygous alt) carrying ``prevalence_percent``
* ``variant combination frequency`` (leaf) -> per population, the joint diplotype
                                       across *both* variants, with an *estimated*
                                       ``prevalence_percent`` (see notes below)

Leaf nodes carry ``prevalence_percent`` so a knowledge-graph query such as
"leaf nodes where prevalence_percent > 20" works directly.

Important modelling notes
-------------------------
* ``rs1800562`` and ``rs1799945`` are *variants* (SNPs) of the HFE gene, not genes.
    - rs1800562 -> HFE p.Cys282Tyr
    - rs1799945 -> HFE p.His63Asp
* "% if alone" is the per-variant ``genotype frequency`` leaf.
* "% for each combo" is the ``variant combination frequency`` leaf. The 1000
  Genomes REST endpoint returns *per-variant* frequencies, not the joint
  two-locus (diplotype) frequencies, so the combination prevalence is
  **estimated under an independence assumption** (product of the two
  per-variant genotype frequencies) and is flagged as such via
  ``estimation_method``. The true joint frequencies require joint genotype
  data and may differ from this estimate.
"""

import json
import logging
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import requests

logger = logging.getLogger(__name__)

ENSEMBL_REST_SERVER = "https://rest.ensembl.org"

# HFE gene metadata (GRCh38, Ensembl).
HFE_GENE = {
    "id": "ENSG00000010704",
    "symbol": "HFE",
    "name": "homeostatic iron regulator",
    "chromosome": "6",
    "start_position": 26087226,
    "end_position": 26098343,
}

# Known amino-acid consequence annotation for the default variants.
VARIANT_CONSEQUENCE = {
    "rs1800562": "HFE p.Cys282Tyr (missense)",
    "rs1799945": "HFE p.His63Asp (missense)",
}


def _sanitize(value: str) -> str:
    """
    AI-Generated.
    Make a string safe to use inside a node identifier.
    """
    return value.replace(":", "_").replace("|", "_").replace("/", "_").replace(" ", "_")


class HfeVariantAdapter:
    """Adapter for HFE variant population prevalence from Ensembl / 1000 Genomes."""

    def __init__(
        self,
        variants: Optional[List[str]] = None,
        gene: Optional[Dict[str, Any]] = None,
        population_prefix: str = "1000GENOMES:phase_3:",
        cache_dir: str | Path = "data",
        raw_data: Optional[Dict[str, Dict[str, Any]]] = None,
        use_cache: bool = True,
        estimate_combinations: bool = True,
        **kwargs: Any,
    ) -> None:
        """
        AI-Generated.

        Args:
        variants: rsIDs to model (default ``["rs1800562", "rs1799945"]``).
        gene: Gene metadata dict (defaults to HFE).
        population_prefix: Only populations whose name starts with this
        prefix are used. Default keeps 1000 Genomes phase 3 populations.
        cache_dir: Directory used to cache Ensembl responses.
        raw_data: Optional pre-fetched ``{rsid: ensembl_response}`` mapping.
        When provided no network calls are made (used for testing).
        use_cache: Whether to read/write cached Ensembl responses on disk.
        estimate_combinations: Whether to emit the joint ``variant
        combination frequency`` leaf nodes.
        **kwargs: Ignored extra configuration.

        """
        self.variants = variants or ["rs1800562", "rs1799945"]
        self.gene = gene or dict(HFE_GENE)
        self.population_prefix = population_prefix
        self.cache_dir = Path(cache_dir)
        self.use_cache = use_cache
        self.estimate_combinations = estimate_combinations
        self.config = kwargs

        self._raw_data: Dict[str, Dict[str, Any]] = raw_data or {}
        # Parsed per-variant structure: rsid -> parsed dict
        self._parsed: Dict[str, Dict[str, Any]] = {}

        logger.info(
            "Initialized HfeVariantAdapter for gene=%s variants=%s",
            self.gene.get("symbol"),
            self.variants,
        )

    # ------------------------------------------------------------------
    # Data acquisition
    # ------------------------------------------------------------------
    def _cache_path(self, rsid: str) -> Path:
        """
        AI-Generated.
        """
        return self.cache_dir / f"ensembl_{rsid}.json"

    def _fetch_variant(self, rsid: str) -> Dict[str, Any]:
        """
        AI-Generated.
        Fetch (or load from cache) the Ensembl variation record for an rsID.
        """
        if rsid in self._raw_data:
            return self._raw_data[rsid]

        cache_path = self._cache_path(rsid)
        if self.use_cache and cache_path.exists():
            logger.info("Loading cached Ensembl data for %s", rsid)
            with open(cache_path, "r", encoding="utf-8") as fh:
                return json.load(fh)

        url = f"{ENSEMBL_REST_SERVER}/variation/human/{rsid}"
        params = {"pops": 1, "genotypes": 1}
        headers = {"Content-Type": "application/json"}
        logger.info("Fetching Ensembl data for %s from %s", rsid, url)
        response = requests.get(url, params=params, headers=headers, timeout=60)
        response.raise_for_status()
        data = response.json()

        if self.use_cache:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            with open(cache_path, "w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=2)

        self._raw_data[rsid] = data
        return data

    # ------------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------------
    def _filter_pop(self, population: str) -> bool:
        """
        AI-Generated.
        """
        return population.startswith(self.population_prefix)

    def _reference_allele(
        self, data: Dict[str, Any], allele_freqs: Dict[str, Dict[str, float]]
    ) -> str:
        """
        AI-Generated.
        Determine the reference allele for zygosity labelling.

        Uses the ancestral allele when available, otherwise the globally most
        frequent allele across the retained populations.

        """
        ancestral = data.get("ancestral_allele")
        observed_alleles = {
            allele for freqs in allele_freqs.values() for allele in freqs
        }
        if ancestral and ancestral in observed_alleles:
            return ancestral

        # Fall back to the most frequent allele overall.
        totals: Dict[str, float] = {}
        for freqs in allele_freqs.values():
            for allele, freq in freqs.items():
                totals[allele] = totals.get(allele, 0.0) + freq
        if totals:
            return max(totals, key=totals.get)
        return ""

    @staticmethod
    def _split_genotype(genotype: str) -> List[str]:
        """
        AI-Generated.
        """
        for sep in ("|", "/"):
            if sep in genotype:
                return genotype.split(sep)
        # Unphased concatenation like "AA".
        return list(genotype)

    def _zygosity(self, alleles: List[str], reference: str) -> str:
        """
        AI-Generated.
        """
        if len(set(alleles)) == 1:
            return (
                "homozygous_reference"
                if alleles[0] == reference
                else "homozygous_alternate"
            )
        return "heterozygous"

    @staticmethod
    def _canonical_genotype(alleles: List[str], reference: str) -> str:
        """
        AI-Generated.
        Return a canonical genotype string (reference allele first).

        Ensures phased/unphased order variants collapse, e.g. ``T|C`` and
        ``C|T`` both become ``C|T`` when ``C`` is the reference allele.

        """
        ordered = sorted(alleles, key=lambda a: (a != reference, a))
        return "|".join(ordered)

    def _observed_genotype_counts(
        self,
        data: Dict[str, Any],
    ) -> Tuple[Dict[Tuple[str, str], int], Dict[str, int], int]:
        """
        AI-Generated.
        Tally per-sample genotypes for the configured population cohort.

        Ensembl returns individual sample genotypes (``genotypes`` field). We
        aggregate the samples belonging to the configured population prefix
        (e.g. all 1000 Genomes phase 3 samples).

        Returns a tuple of:
        * genotype tally keyed by a sorted allele pair,
        * per-allele tally (to determine reference/alternate empirically),
        * total number of samples counted.

        """
        genotype_counts: Dict[Tuple[str, str], int] = {}
        allele_counts: Dict[str, int] = {}
        total = 0
        for entry in data.get("genotypes", []) or []:
            sample = entry.get("sample", "")
            if not sample or not sample.startswith(self.population_prefix):
                continue
            genotype = entry.get("genotype")
            if not genotype:
                continue
            alleles = self._split_genotype(genotype)
            if len(alleles) != 2 or any(a in (".", "-", "") for a in alleles):
                continue
            key = tuple(sorted(alleles))
            genotype_counts[key] = genotype_counts.get(key, 0) + 1
            for a in alleles:
                allele_counts[a] = allele_counts.get(a, 0) + 1
            total += 1
        return genotype_counts, allele_counts, total

    def _build_observed_frequencies(
        self,
        genotype_counts: Dict[Tuple[str, str], int],
        total: int,
        reference: str,
    ) -> Dict[str, Dict[str, Any]]:
        """
        AI-Generated.
        Turn observed genotype counts into a genotype-frequency table.

        Handles multi-allelic sites by classifying every observed genotype
        rather than assuming a single alternate allele.

        """
        genotypes: Dict[str, Dict[str, Any]] = {}
        for pair, count in genotype_counts.items():
            alleles = self._canonical_genotype(list(pair), reference).split("|")
            canon = "|".join(alleles)
            genotypes[canon] = {
                "genotype": canon,
                "alleles": alleles,
                "zygosity": self._zygosity(alleles, reference),
                "frequency": round(count / total, 6),
                "count": count,
                "sample_size": total,
                "method": "observed_genotype_counts",
            }
        return genotypes

    def _parse_variant(self, rsid: str, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        AI-Generated.
        Parse an Ensembl response into allele and genotype frequency tables.
        """
        # ---- allele frequencies per (retained) population -------------
        allele_freqs: Dict[str, Dict[str, float]] = {}
        allele_counts: Dict[str, Dict[str, int]] = {}
        for entry in data.get("populations", []):
            pop = entry.get("population")
            if not pop or not self._filter_pop(pop):
                continue
            allele = entry.get("allele")
            freq = entry.get("frequency")
            if allele is None or freq is None:
                continue
            allele_freqs.setdefault(pop, {})[allele] = float(freq)
            if entry.get("allele_count") is not None:
                allele_counts.setdefault(pop, {})[allele] = int(entry["allele_count"])

        # ---- observed genotype tally (drives reference/alt inference) -
        mapping = (data.get("mappings") or [{}])[0]
        allele_string = mapping.get("allele_string") or ""
        allele_string_parts = [a for a in allele_string.split("/") if a]

        genotype_counts, obs_allele_counts, obs_total = self._observed_genotype_counts(
            data
        )

        # Reference allele: the genomic reference (allele_string[0]) when it is
        # actually observed, otherwise the most frequent observed allele.
        if allele_string_parts and allele_string_parts[0] in obs_allele_counts:
            reference = allele_string_parts[0]
        elif obs_allele_counts:
            reference = max(obs_allele_counts, key=obs_allele_counts.get)
        elif allele_string_parts:
            reference = allele_string_parts[0]
        else:
            reference = self._reference_allele(data, allele_freqs)

        # Alternate allele: most frequent observed non-reference allele.
        non_ref = {a: c for a, c in obs_allele_counts.items() if a != reference}
        if non_ref:
            alt = max(non_ref, key=non_ref.get)
        else:
            alts = sorted(
                {a for f in allele_freqs.values() for a in f if a != reference}
            )
            alt = alts[0] if alts else None

        genotype_freqs: Dict[str, Dict[str, Dict[str, Any]]] = {}

        # ---- observed aggregate genotype frequencies (real counts) ----
        observed_pop = f"{self.population_prefix}ALL"
        if obs_total > 0:
            genotype_freqs[observed_pop] = self._build_observed_frequencies(
                genotype_counts, obs_total, reference
            )

        # ---- per-population Hardy-Weinberg estimates ------------------
        for pop, freqs in allele_freqs.items():
            if pop == observed_pop:
                # Prefer the real observed counts over an HWE estimate.
                continue
            if reference not in freqs:
                continue
            alts = [a for a in freqs if a != reference]
            if len(alts) != 1:
                # Skip non-biallelic populations for the HWE estimate.
                continue
            alt_allele = alts[0]
            p = freqs[reference]
            q = freqs[alt_allele]
            genotype_freqs[pop] = {
                f"{reference}|{reference}": {
                    "genotype": f"{reference}|{reference}",
                    "alleles": [reference, reference],
                    "zygosity": "homozygous_reference",
                    "frequency": round(p * p, 6),
                    "method": "hardy_weinberg",
                },
                f"{reference}|{alt_allele}": {
                    "genotype": f"{reference}|{alt_allele}",
                    "alleles": [reference, alt_allele],
                    "zygosity": "heterozygous",
                    "frequency": round(2 * p * q, 6),
                    "method": "hardy_weinberg",
                },
                f"{alt_allele}|{alt_allele}": {
                    "genotype": f"{alt_allele}|{alt_allele}",
                    "alleles": [alt_allele, alt_allele],
                    "zygosity": "homozygous_alternate",
                    "frequency": round(q * q, 6),
                    "method": "hardy_weinberg",
                },
            }

        parsed = {
            "rsid": rsid,
            "reference_allele": reference,
            "alternate_allele": alt,
            "minor_allele": data.get("minor_allele"),
            "ancestral_allele": mapping.get("ancestral_allele")
            or data.get("ancestral_allele"),
            "maf": data.get("MAF"),
            "consequence": VARIANT_CONSEQUENCE.get(rsid),
            "allele_string": allele_string or None,
            "chromosome": mapping.get("seq_region_name"),
            "position": mapping.get("start"),
            "assembly": mapping.get("assembly_name"),
            "allele_freqs": allele_freqs,
            "allele_counts": allele_counts,
            "genotype_freqs": genotype_freqs,
        }
        return parsed

    def _ensure_parsed(self) -> None:
        """
        AI-Generated.
        """
        if self._parsed:
            return
        for rsid in self.variants:
            data = self._fetch_variant(rsid)
            self._parsed[rsid] = self._parse_variant(rsid, data)

    # ------------------------------------------------------------------
    # Node / edge id helpers
    # ------------------------------------------------------------------
    def _gene_id(self) -> str:
        """
        AI-Generated.
        """
        return self.gene["id"]

    @staticmethod
    def _population_id(pop: str) -> str:
        """
        AI-Generated.
        """
        return _sanitize(pop)

    @staticmethod
    def _population_label(pop: str) -> str:
        # e.g. "1000GENOMES:phase_3:GBR" -> "GBR"
        """
        AI-Generated.
        """
        return pop.split(":")[-1]

    def _genotype_freq_id(self, rsid: str, pop: str, genotype: str) -> str:
        """
        AI-Generated.
        """
        return f"gtf_{rsid}_{_sanitize(pop)}_{_sanitize(genotype)}"

    def _combination_id(self, pop: str, pairs: List[Tuple[str, str]]) -> str:
        """
        AI-Generated.
        """
        combo = "__".join(f"{rsid}_{_sanitize(gt)}" for rsid, gt in pairs)
        return f"combo_{_sanitize(pop)}_{combo}"

    # ------------------------------------------------------------------
    # BioCypher interface
    # ------------------------------------------------------------------
    def get_nodes(self) -> Iterable[Tuple[str, str, Dict[str, Any]]]:
        """
        AI-Generated.
        Yield ``(node_id, node_label, properties)`` tuples.
        """
        self._ensure_parsed()
        seen_populations: set[str] = set()

        # ---- gene node -----------------------------------------------
        yield (
            self._gene_id(),
            "gene",
            {
                "name": self.gene.get("name"),
                "symbol": self.gene.get("symbol"),
                "chromosome": self.gene.get("chromosome"),
                "start_position": self.gene.get("start_position"),
                "end_position": self.gene.get("end_position"),
            },
        )

        for rsid, parsed in self._parsed.items():
            # ---- variant node ----------------------------------------
            yield (
                rsid,
                "sequence variant",
                {
                    "rsid": rsid,
                    "consequence": parsed.get("consequence"),
                    "reference_allele": parsed.get("reference_allele"),
                    "alternate_allele": parsed.get("alternate_allele"),
                    "ancestral_allele": parsed.get("ancestral_allele"),
                    "minor_allele": parsed.get("minor_allele"),
                    "maf": parsed.get("maf"),
                    "allele_string": parsed.get("allele_string"),
                    "chromosome": parsed.get("chromosome"),
                    "position": parsed.get("position"),
                    "assembly": parsed.get("assembly"),
                },
            )

            # ---- population + genotype frequency leaf nodes ----------
            for pop, genotypes in parsed["genotype_freqs"].items():
                if pop not in seen_populations:
                    seen_populations.add(pop)
                    label = self._population_label(pop)
                    yield (
                        self._population_id(pop),
                        "population of individual organisms",
                        {
                            "name": label,
                            "population_code": pop,
                            "is_aggregate": label
                            in {"ALL", "AFR", "AMR", "EAS", "EUR", "SAS"},
                            "source": "1000 Genomes Project",
                        },
                    )

                for genotype, info in genotypes.items():
                    prevalence = round(info["frequency"] * 100, 4)
                    props = {
                        "rsid": rsid,
                        "population_code": pop,
                        "genotype": genotype,
                        "zygosity": info["zygosity"],
                        "alleles": info["alleles"],
                        # "prevalence if alone" for this single variant
                        "prevalence_percent": prevalence,
                        "frequency": info["frequency"],
                        "estimation_method": info.get("method", "unknown"),
                    }
                    if "count" in info:
                        props["genotype_count"] = info["count"]
                        props["sample_size"] = info["sample_size"]
                    yield (
                        self._genotype_freq_id(rsid, pop, genotype),
                        "genotype frequency",
                        props,
                    )

        # ---- variant combination (joint diplotype) leaf nodes --------
        if self.estimate_combinations:
            yield from self._combination_nodes()

    def _combination_nodes(
        self,
    ) -> Iterable[Tuple[str, str, Dict[str, Any]]]:
        """
        AI-Generated.
        Yield joint-diplotype leaf nodes estimated under independence.
        """
        rsids = list(self._parsed.keys())
        if len(rsids) < 2:
            return

        # Populations common to all variants.
        common_pops = None
        for rsid in rsids:
            pops = set(self._parsed[rsid]["genotype_freqs"].keys())
            common_pops = pops if common_pops is None else common_pops & pops
        common_pops = common_pops or set()

        for pop in sorted(common_pops):
            # Build the cartesian product of genotypes across all variants.
            per_variant = [
                [
                    (rsid, gt, info)
                    for gt, info in self._parsed[rsid]["genotype_freqs"][pop].items()
                ]
                for rsid in rsids
            ]

            for combo in _product(per_variant):
                pairs = [(rsid, gt) for rsid, gt, _ in combo]
                joint_freq = 1.0
                for _, _, info in combo:
                    joint_freq *= info["frequency"]
                prevalence = round(joint_freq * 100, 4)

                props: Dict[str, Any] = {
                    "population_code": pop,
                    "prevalence_percent": prevalence,
                    "frequency": round(joint_freq, 8),
                    "estimation_method": "independence_assumption",
                    "note": (
                        "Joint diplotype prevalence estimated as the product of "
                        "per-variant genotype frequencies under independence; "
                        "the true joint value may differ."
                    ),
                }
                for rsid, gt, info in combo:
                    props[f"{rsid}_genotype"] = gt
                    props[f"{rsid}_zygosity"] = info["zygosity"]
                combo_label = " + ".join(f"{rsid}:{gt}" for rsid, gt, _ in combo)
                props["combination"] = combo_label

                yield (
                    self._combination_id(pop, pairs),
                    "variant combination frequency",
                    props,
                )

    def get_edges(
        self,
    ) -> Iterable[Tuple[str, str, str, str, Dict[str, Any]]]:
        """
        AI-Generated.
        Yield ``(edge_id, source_id, target_id, label, properties)`` tuples.
        """
        self._ensure_parsed()
        gene_id = self._gene_id()

        for rsid, parsed in self._parsed.items():
            # gene -> variant
            yield (
                f"has_variant_{gene_id}_{rsid}",
                gene_id,
                rsid,
                "has_variant",
                {},
            )

            for pop, genotypes in parsed["genotype_freqs"].items():
                pop_id = self._population_id(pop)
                for genotype in genotypes:
                    gtf_id = self._genotype_freq_id(rsid, pop, genotype)
                    # variant -> genotype frequency
                    yield (
                        f"has_gtf_{gtf_id}",
                        rsid,
                        gtf_id,
                        "has_genotype_frequency",
                        {},
                    )
                    # genotype frequency -> population
                    yield (
                        f"gtf_in_pop_{gtf_id}",
                        gtf_id,
                        pop_id,
                        "genotype_frequency_in_population",
                        {},
                    )

        if self.estimate_combinations:
            yield from self._combination_edges()

    def _combination_edges(
        self,
    ) -> Iterable[Tuple[str, str, str, str, Dict[str, Any]]]:
        """
        AI-Generated.
        """
        rsids = list(self._parsed.keys())
        if len(rsids) < 2:
            return

        common_pops = None
        for rsid in rsids:
            pops = set(self._parsed[rsid]["genotype_freqs"].keys())
            common_pops = pops if common_pops is None else common_pops & pops
        common_pops = common_pops or set()

        for pop in sorted(common_pops):
            pop_id = self._population_id(pop)
            per_variant = [
                [
                    (rsid, gt, info)
                    for gt, info in self._parsed[rsid]["genotype_freqs"][pop].items()
                ]
                for rsid in rsids
            ]
            for combo in _product(per_variant):
                pairs = [(rsid, gt) for rsid, gt, _ in combo]
                combo_id = self._combination_id(pop, pairs)
                # combination -> population
                yield (
                    f"combo_in_pop_{combo_id}",
                    combo_id,
                    pop_id,
                    "combination_in_population",
                    {},
                )
                # combination -> each involved variant
                for rsid, _gt in pairs:
                    yield (
                        f"combo_involves_{combo_id}_{rsid}",
                        combo_id,
                        rsid,
                        "combination_involves_variant",
                        {},
                    )

    # ------------------------------------------------------------------
    # Utility
    # ------------------------------------------------------------------
    def get_metadata(self) -> Dict[str, Any]:
        """
        AI-Generated.
        """
        return {
            "name": "HfeVariantAdapter",
            "gene": self.gene.get("symbol"),
            "variants": self.variants,
            "data_source": "Ensembl REST API (1000 Genomes Project)",
            "data_type": "api",
            "version": "0.1.0",
            "adapter_class": "HfeVariantAdapter",
        }

    def validate_data_source(self) -> bool:
        """
        AI-Generated.
        Validate that variant data can be parsed for all requested variants.
        """
        try:
            self._ensure_parsed()
            return all(
                bool(self._parsed[rsid]["genotype_freqs"]) for rsid in self.variants
            )
        except Exception as exc:  # pragma: no cover - defensive
            logger.error("Data source validation failed: %s", exc)
            return False


def _product(list_of_lists: List[List[Any]]) -> Iterable[Tuple[Any, ...]]:
    """
    AI-Generated.
    Cartesian product helper (itertools.product wrapper for readability).
    """
    from itertools import product

    yield from product(*list_of_lists)
