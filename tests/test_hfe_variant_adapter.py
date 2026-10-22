"""
Tests for HfeVariantAdapter.

These tests use injected ``raw_data`` (mimicking real Ensembl REST responses,
including per-sample genotypes) so no network access is required and results are
deterministic.
"""

import math

import pytest

from hfe_variant_prevalence.adapters.hfe_variant_adapter import HfeVariantAdapter

PREFIX = "1000GENOMES:phase_3:"


def _ensembl_like(rsid, ref, alt, sample_genotypes, pop_allele_freqs):
    """
    AI-Generated.
    Build a minimal Ensembl-like response.

    sample_genotypes: list of genotype strings (one per 1000G sample), e.g.
    ["C|C", "C|T", "T|T", ...]
    pop_allele_freqs: {population: {allele: frequency}}  (sub/super populations)

    """
    genotypes = []
    for i, gt in enumerate(sample_genotypes):
        genotypes.append(
            {"genotype": gt, "sample": f"{PREFIX}HG{i:05d}", "gender": "Male"}
        )
    # Add a non-1000G sample that must be ignored by the aggregate count.
    genotypes.append({"genotype": "T|T", "sample": "OTHER:sampleX"})

    populations = []
    for pop, freqs in pop_allele_freqs.items():
        for allele, freq in freqs.items():
            populations.append({"population": pop, "allele": allele, "frequency": freq})
    # Add a gnomAD population that must be filtered out.
    populations.append({"population": "gnomADg:ALL", "allele": ref, "frequency": 0.9})
    populations.append({"population": "gnomADg:ALL", "allele": alt, "frequency": 0.1})

    return {
        "name": rsid,
        "minor_allele": alt,
        "MAF": 0.4,
        "mappings": [
            {
                "seq_region_name": "6",
                "start": 26092913,
                "allele_string": f"{ref}/{alt}",
                "assembly_name": "GRCh38",
                "ancestral_allele": alt,
            }
        ],
        "populations": populations,
        "genotypes": genotypes,
    }


@pytest.fixture
def raw_data():
    """
    AI-Generated.
    """
    return {
        # rs1800562: 10 samples -> 4 C|C, 4 C|T, 2 T|T  (mix phased order)
        "rs1800562": _ensembl_like(
            "rs1800562",
            "C",
            "T",
            ["C|C", "C|C", "C|C", "C|C", "C|T", "T|C", "C|T", "C|T", "T|T", "T|T"],
            {
                f"{PREFIX}GBR": {"C": 0.6, "T": 0.4},
                f"{PREFIX}EUR": {"C": 0.7, "T": 0.3},
            },
        ),
        # rs1799945: 10 samples -> 2 C|C, 5 C|T, 3 T|T
        "rs1799945": _ensembl_like(
            "rs1799945",
            "C",
            "T",
            ["C|C", "C|C", "C|T", "C|T", "C|T", "C|T", "C|T", "T|T", "T|T", "T|T"],
            {
                f"{PREFIX}GBR": {"C": 0.5, "T": 0.5},
                f"{PREFIX}EUR": {"C": 0.55, "T": 0.45},
            },
        ),
    }


@pytest.fixture
def adapter(raw_data):
    """
    AI-Generated.
    """
    return HfeVariantAdapter(raw_data=raw_data, use_cache=False)


class TestHfeVariantAdapter:
    def test_metadata(self, adapter):
        """
        AI-Generated.
        """
        md = adapter.get_metadata()
        assert md["adapter_class"] == "HfeVariantAdapter"
        assert md["gene"] == "HFE"
        assert md["variants"] == ["rs1800562", "rs1799945"]
        assert md["data_type"] == "api"

    def test_validate_data_source(self, adapter):
        """
        AI-Generated.
        """
        assert adapter.validate_data_source() is True

    def test_node_tuple_shape(self, adapter):
        """
        AI-Generated.
        """
        nodes = list(adapter.get_nodes())
        assert len(nodes) > 0
        for node_id, label, props in nodes:
            assert isinstance(node_id, str)
            assert isinstance(label, str)
            assert isinstance(props, dict)

    def test_edge_tuple_shape(self, adapter):
        """
        AI-Generated.
        """
        edges = list(adapter.get_edges())
        assert len(edges) > 0
        for edge in edges:
            assert len(edge) == 5
            eid, src, tgt, label, props = edge
            assert all(isinstance(x, str) for x in (eid, src, tgt, label))
            assert isinstance(props, dict)

    def test_all_expected_node_labels(self, adapter):
        """
        AI-Generated.
        """
        labels = {n[1] for n in adapter.get_nodes()}
        assert labels == {
            "gene",
            "sequence variant",
            "population of individual organisms",
            "genotype frequency",
            "variant combination frequency",
        }

    def test_variant_nodes(self, adapter):
        """
        AI-Generated.
        """
        nodes = list(adapter.get_nodes())
        variants = {n[0]: n[2] for n in nodes if n[1] == "sequence variant"}
        assert set(variants) == {"rs1800562", "rs1799945"}
        assert variants["rs1800562"]["reference_allele"] == "C"
        assert variants["rs1800562"]["alternate_allele"] == "T"
        assert "Cys282Tyr" in variants["rs1800562"]["consequence"]

    def test_population_filtering(self, adapter):
        """
        AI-Generated.
        """
        pops = {
            n[2]["population_code"]
            for n in adapter.get_nodes()
            if n[1] == "population of individual organisms"
        }
        # gnomAD and OTHER filtered out; aggregate ALL added from samples.
        assert pops == {
            f"{PREFIX}ALL",
            f"{PREFIX}GBR",
            f"{PREFIX}EUR",
        }

    def test_observed_genotype_prevalence(self, adapter):
        """
        AI-Generated.
        Real observed counts for the aggregate 1000G cohort.
        """
        gtf = [
            n[2]
            for n in adapter.get_nodes()
            if n[1] == "genotype frequency"
            and n[2]["population_code"] == f"{PREFIX}ALL"
            and n[2]["rsid"] == "rs1800562"
        ]
        by_zyg = {g["zygosity"]: g for g in gtf}
        # 4/10 C|C, 4/10 het, 2/10 T|T (non-1000G sample ignored)
        assert by_zyg["homozygous_reference"]["prevalence_percent"] == pytest.approx(
            40.0
        )
        assert by_zyg["heterozygous"]["prevalence_percent"] == pytest.approx(40.0)
        assert by_zyg["homozygous_alternate"]["prevalence_percent"] == pytest.approx(
            20.0
        )
        assert all(g["estimation_method"] == "observed_genotype_counts" for g in gtf)
        assert by_zyg["heterozygous"]["sample_size"] == 10
        assert by_zyg["heterozygous"]["genotype_count"] == 4
        # Phased order T|C collapsed into canonical C|T
        assert by_zyg["heterozygous"]["genotype"] == "C|T"

    def test_observed_frequencies_sum_to_one(self, adapter):
        """
        AI-Generated.
        """
        for rsid in ("rs1800562", "rs1799945"):
            freqs = [
                n[2]["frequency"]
                for n in adapter.get_nodes()
                if n[1] == "genotype frequency"
                and n[2]["population_code"] == f"{PREFIX}ALL"
                and n[2]["rsid"] == rsid
            ]
            assert math.isclose(sum(freqs), 1.0, rel_tol=1e-6)

    def test_hardy_weinberg_subpopulations(self, adapter):
        """
        AI-Generated.
        """
        gtf = [
            n[2]
            for n in adapter.get_nodes()
            if n[1] == "genotype frequency"
            and n[2]["population_code"] == f"{PREFIX}EUR"
            and n[2]["rsid"] == "rs1800562"
        ]
        by_zyg = {g["zygosity"]: g for g in gtf}
        # p=0.7 -> p^2 = 0.49 -> 49%
        assert by_zyg["homozygous_reference"]["prevalence_percent"] == pytest.approx(
            49.0
        )
        assert by_zyg["heterozygous"]["prevalence_percent"] == pytest.approx(42.0)
        assert by_zyg["homozygous_alternate"]["prevalence_percent"] == pytest.approx(
            9.0
        )
        assert all(g["estimation_method"] == "hardy_weinberg" for g in gtf)

    def test_combination_prevalence_independence(self, adapter):
        """
        AI-Generated.
        """
        combos = [
            n[2]
            for n in adapter.get_nodes()
            if n[1] == "variant combination frequency"
            and n[2]["population_code"] == f"{PREFIX}ALL"
        ]
        # 3 genotypes x 3 genotypes = 9 combos for the ALL population
        assert len(combos) == 9
        # rs1800562 T|T (0.2) x rs1799945 T|T (0.3) = 0.06 -> 6%
        target = next(
            c
            for c in combos
            if c["rs1800562_genotype"] == "T|T" and c["rs1799945_genotype"] == "T|T"
        )
        assert target["prevalence_percent"] == pytest.approx(6.0)
        assert target["estimation_method"] == "independence_assumption"

    def test_combination_probabilities_sum_to_one(self, adapter):
        """
        AI-Generated.
        """
        freqs = [
            n[2]["frequency"]
            for n in adapter.get_nodes()
            if n[1] == "variant combination frequency"
            and n[2]["population_code"] == f"{PREFIX}ALL"
        ]
        assert math.isclose(sum(freqs), 1.0, rel_tol=1e-6)

    def test_edge_labels(self, adapter):
        """
        AI-Generated.
        """
        labels = {e[3] for e in adapter.get_edges()}
        assert labels == {
            "has_variant",
            "has_genotype_frequency",
            "genotype_frequency_in_population",
            "combination_in_population",
            "combination_involves_variant",
        }
        gene_edges = [e for e in adapter.get_edges() if e[3] == "has_variant"]
        assert {e[2] for e in gene_edges} == {"rs1800562", "rs1799945"}

    def test_leaf_prevalence_query(self, adapter):
        """
        AI-Generated.
        Emulate the target use-case: leaf nodes with prevalence > 20%.
        """
        leaves = [
            (n[1], n[2])
            for n in adapter.get_nodes()
            if n[1] in ("genotype frequency", "variant combination frequency")
            and n[2]["prevalence_percent"] > 20
        ]
        assert leaves, "expected at least one leaf above 20% prevalence"
        assert all(props["prevalence_percent"] > 20 for _, props in leaves)
