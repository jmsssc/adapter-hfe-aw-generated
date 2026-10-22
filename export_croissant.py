"""Export the bundled prevalence data as CSV tables and Croissant metadata."""

import csv
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

import mlcroissant as mlc

from hfe_variant_prevalence.adapters.hfe_variant_adapter import HfeVariantAdapter

ROOT = Path(__file__).resolve().parent
VARIANTS = ("rs1800562", "rs1799945")


def write_table(
    name: str,
    rows: Sequence[Mapping[str, object]],
    columns: Mapping[str, str],
    description: str,
) -> dict[str, object]:
    """
    AI-Generated.
    Write one CSV table and describe its fields for Croissant readers.
    """
    path = ROOT / "data" / f"{name}.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    file_id = path.relative_to(ROOT).as_posix()
    return {
        "@type": "cr:RecordSet",
        "@id": name,
        "name": name,
        "description": description,
        "key": {"@id": f"{name}/id"},
        "field": [
            {
                "@type": "cr:Field",
                "@id": f"{name}/{column}",
                "name": column,
                "description": column.replace("_", " "),
                "dataType": data_type,
                "source": {
                    "fileObject": {"@id": file_id},
                    "extract": {"column": column},
                },
            }
            for column, data_type in columns.items()
        ],
    }


def main() -> None:
    """
    AI-Generated.
    Describe cached Ensembl data and export prevalence without network requests.
    """
    raw_data = {
        rsid: json.loads((ROOT / "data" / f"ensembl_{rsid}.json").read_text())
        for rsid in VARIANTS
    }
    adapter = HfeVariantAdapter(raw_data=raw_data, use_cache=False)
    nodes = list(adapter.get_nodes())
    genotype_columns = {
        "id": "sc:Text",
        "rsid": "sc:Text",
        "population_code": "sc:Text",
        "genotype": "sc:Text",
        "zygosity": "sc:Text",
        "prevalence_percent": "sc:Float",
        "frequency": "sc:Float",
        "estimation_method": "sc:Text",
    }
    combination_columns = {
        "id": "sc:Text",
        "population_code": "sc:Text",
        "combination": "sc:Text",
        "prevalence_percent": "sc:Float",
        "frequency": "sc:Float",
        "estimation_method": "sc:Text",
        "note": "sc:Text",
    }
    record_sets = []
    for label, name, columns, description in (
        (
            "genotype frequency",
            "genotype_frequencies",
            genotype_columns,
            "Per-variant genotype prevalence. Aggregate ALL values use observed "
            "sample counts; other populations use Hardy-Weinberg estimates.",
        ),
        (
            "variant combination frequency",
            "variant_combinations",
            combination_columns,
            "Joint diplotype estimates calculated by multiplying per-variant "
            "frequencies under independence; these are not observed joint counts.",
        ),
    ):
        rows = [
            {"id": node_id, **properties}
            for node_id, node_label, properties in nodes
            if node_label == label
        ]
        record_sets.append(write_table(name, rows, columns, description))

    distribution = []
    paths = sorted((ROOT / "data").glob("*"))
    paths.extend(sorted((ROOT / "biocypher-out").glob("*.csv")))
    for path in paths:
        if not path.is_file() or path.suffix not in {".json", ".csv"}:
            continue
        relative = path.relative_to(ROOT).as_posix()
        encoding = "application/json" if path.suffix == ".json" else "text/csv"
        distribution.append(
            {
                "@type": "cr:FileObject",
                "@id": relative,
                "name": path.name,
                "contentUrl": relative,
                "encodingFormat": encoding,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    metadata = {
        "@context": mlc.Metadata(name="context").to_json()["@context"],
        "@type": "sc:Dataset",
        "name": "hfe-variant-prevalence",
        "description": "HFE rs1800562/rs1799945 population prevalence from cached "
        "Ensembl REST responses for the 1000 Genomes Project. Includes raw "
        "responses, BioCypher graph CSVs, and Croissant-readable prevalence "
        "tables. Subpopulation genotypes use Hardy-Weinberg estimates. Joint "
        "diplotypes assume independence and are not observed joint prevalence.",
        "conformsTo": "http://mlcommons.org/croissant/1.0",
        "version": "0.1.0",
        "keywords": ["HFE", "1000 Genomes", "population prevalence", "BioCypher"],
        "distribution": distribution,
        "recordSet": record_sets,
    }
    (ROOT / "croissant.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Exported {len(record_sets)} record sets to {ROOT / 'croissant.json'}")


if __name__ == "__main__":
    main()
