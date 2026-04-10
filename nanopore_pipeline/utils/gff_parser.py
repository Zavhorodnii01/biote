"""
GFF3 Ground Truth Parser
=========================
Parses a GFF3 annotation file (from NCBI RefSeq) to extract known gene names
and match them against VFDB/CARD category keyword lists.

Used in Test 4 (real reference genome) to build a ground-truth gene set so the
pipeline can report "found X / Y annotated genes" per category.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from nanopore_pipeline.utils.logging_config import setup_logging

logger = setup_logging("utils.gff_parser")


class GFFGroundTruth:
    """
    Parses GFF3 annotation and matches genes to pathogenicity categories.

    Usage:
        gt = GFFGroundTruth()
        ground_truth = gt.build(gff_path, vf_categories, amr_categories)
        # ground_truth = {
        #   "Exotoxin": ["hla", "lukF", ...],
        #   "Beta-Lactam Resistance": ["mecA", "blaZ", ...],
        #   ...
        # }
    """

    # GFF9 attribute keys that may contain a gene name / product description
    _NAME_KEYS = ("gene", "Name", "product", "gene_synonym", "locus_tag", "note")

    def build(
        self,
        gff_path: Path,
        vf_categories: dict[str, list[str]],
        amr_categories: dict[str, list[str]],
    ) -> dict[str, list[str]]:
        """
        Parse the GFF3 file and return a dict mapping category → list of gene
        names found in the annotation that match that category.
        """
        genes = self.parse_gff(gff_path)
        logger.info("GFF parsed: %d gene/product entries found in %s", len(genes), gff_path.name)

        all_categories: dict[str, list[str]] = {}
        all_categories.update(self.match_to_categories(genes, vf_categories))
        all_categories.update(self.match_to_categories(genes, amr_categories))

        total = sum(len(v) for v in all_categories.values())
        logger.info("Ground truth: %d genes mapped across %d categories", total, len(all_categories))
        return all_categories

    def parse_gff(self, gff_path: Path) -> list[dict]:
        """
        Read a GFF3 file and extract gene-related records.
        Returns a list of dicts with keys: seqid, type, start, end, attributes_raw,
        gene_name, product.
        """
        genes: list[dict] = []
        path = Path(gff_path)
        if not path.exists():
            logger.warning("GFF file not found: %s", path)
            return genes

        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue

                feature_type = parts[2].lower()
                # Focus on gene / CDS / mRNA records (skip region/repeat etc.)
                if feature_type not in ("gene", "cds", "mrna", "rna"):
                    continue

                attrs = self._parse_attributes(parts[8])
                gene_name = attrs.get("gene") or attrs.get("Name") or attrs.get("locus_tag", "")
                product = attrs.get("product", "")

                if not gene_name and not product:
                    continue

                genes.append({
                    "seqid": parts[0],
                    "type": feature_type,
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "gene_name": gene_name.lower(),
                    "product": product.lower(),
                    "attrs": attrs,
                })

        return genes

    def match_to_categories(
        self,
        genes: list[dict],
        categories: dict[str, list[str]],
    ) -> dict[str, list[str]]:
        """
        For each category, find which genes match any of the category keywords.
        Returns {category: [matched_gene_names]}.
        """
        result: dict[str, list[str]] = {}
        for category, keywords in categories.items():
            matched: list[str] = []
            seen: set[str] = set()
            for gene in genes:
                searchable = gene["gene_name"] + " " + gene["product"]
                for kw in keywords:
                    if kw.lower() in searchable:
                        label = gene["gene_name"] or gene["product"][:30]
                        if label and label not in seen:
                            matched.append(label)
                            seen.add(label)
                        break
            if matched:
                result[category] = matched
        return result

    @staticmethod
    def _parse_attributes(attr_string: str) -> dict[str, str]:
        """Parse GFF3 column 9 attributes into a key→value dict."""
        result: dict[str, str] = {}
        for part in attr_string.split(";"):
            part = part.strip()
            if "=" in part:
                key, _, value = part.partition("=")
                result[key.strip()] = value.strip().replace("%2C", ",").replace("%3B", ";")
        return result
