"""
Central configuration for the Nanopore Pathogenicity Diagnostic Pipeline.
All paths, thresholds, and external tool settings live here.
"""

import os
from pathlib import Path

# ── Project root ──────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"

# ── Database ──────────────────────────────────────────────────────────────────
DATABASE_URL = os.getenv(
    "NANOPORE_DB_URL",
    f"sqlite:///{BASE_DIR / 'data' / 'nanopore_pipeline.db'}"
)

# ── Reference data paths ─────────────────────────────────────────────────────
VFDB_DIR = DATA_DIR / "vfdb"
CARD_DIR = DATA_DIR / "card"
SAMPLES_DIR = DATA_DIR / "samples"
RESULTS_DIR = DATA_DIR / "results"
ASSEMBLIES_DIR = DATA_DIR / "assemblies"

VFDB_FASTA_URL = "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz"
VFDB_PROTEIN_URL = "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz"
CARD_DATA_URL = "https://card.mcmaster.ca/latest/data"

# ── Alignment settings ───────────────────────────────────────────────────────
BLAST_EVALUE = 1e-10
BLAST_MIN_IDENTITY = 70.0       # percent
BLAST_MIN_COVERAGE = 50.0       # percent
BLAST_NUM_THREADS = os.cpu_count() or 4
BLAST_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

DIAMOND_EVALUE = 1e-10
DIAMOND_MIN_IDENTITY = 70.0
DIAMOND_MIN_COVERAGE = 50.0
DIAMOND_BLOCK_SIZE = 2.0
DIAMOND_INDEX_CHUNKS = 4

# ── Nanopore-specific parameters ─────────────────────────────────────────────
NANOPORE_MIN_READ_LENGTH = 200
# Word size tuned for ~10% Nanopore error rate
BLAST_WORD_SIZE = 7
BLAST_REWARD = 1
BLAST_PENALTY = -1
BLAST_GAPOPEN = 2
BLAST_GAPEXTEND = 1

# ── Assembly settings ─────────────────────────────────────────────────────────
# Flye assembler parameters (tuned for Nanopore metagenomics)
FLYE_GENOME_SIZE = "5m"                # 5 megabases (typical bacterial genome)
FLYE_MIN_OVERLAP = 1000                # Minimum overlap between reads (bp)
FLYE_ITERATIONS = 1                    # Polishing iterations (balance speed/quality)
FLYE_THREADS = os.cpu_count() or 4
FLYE_META_MODE = True                  # Critical for metagenomics (multiple strains)
FLYE_MIN_READ_LENGTH = 1000            # Filter short reads before assembly
FLYE_TIMEOUT = 3600                    # Assembly timeout in seconds (1 hour)

# ── Pathogenicity classification categories ──────────────────────────────────
VF_CATEGORIES = {
    "Adherence": ["adherence", "adhesin", "pilus", "fimbriae", "attachment", "binding"],
    "Motility": ["flagell", "motility", "chemotaxis", "type iv pil"],
    "Exotoxin": ["exotoxin", "toxin", "hemolysin", "cytolysin", "leukocidin", "enterotoxin"],
    "Endotoxin": ["endotoxin", "lipopolysaccharide", "lps", "lipid a"],
    "Secretion System": [
        "type i secretion", "type ii secretion", "type iii secretion",
        "type iv secretion", "type v secretion", "type vi secretion",
        "t1ss", "t2ss", "t3ss", "t4ss", "t5ss", "t6ss",
        "secretion system",
    ],
    "Iron Uptake": ["siderophore", "iron", "heme", "ferric", "yersiniabactin", "aerobactin"],
    "Immune Evasion": ["capsule", "immune evasion", "serum resistance", "complement"],
    "Biofilm": ["biofilm", "quorum sensing", "autoinducer"],
    "Invasion": ["invasion", "invasin", "internalin"],
    "Regulation": ["regulator", "two-component", "transcriptional activator"],
}

AMR_CATEGORIES = {
    "Beta-Lactam Resistance": ["beta-lactam", "oxa", "tem", "shv", "ctx-m", "ndm", "kpc", "vim", "imp"],
    "Aminoglycoside Resistance": ["aminoglycoside", "aac", "ant", "aph"],
    "Fluoroquinolone Resistance": ["fluoroquinolone", "qnr", "gyrase", "topoisomerase"],
    "Tetracycline Resistance": ["tetracycline", "tet("],
    "Macrolide Resistance": ["macrolide", "erm", "mef"],
    "Glycopeptide Resistance": ["vancomycin", "glycopeptide", "van"],
    "Sulfonamide Resistance": ["sulfonamide", "sul1", "sul2"],
    "Efflux Pump": ["efflux", "mex", "acr", "tolc"],
}

# ── Category difficulty (bioinformatic detectability) ────────────────────────
# Based on gene size, copy number, and mutation-level changes.
# Easy   = large genes or multi-gene clusters (hard to miss)
# Medium = moderately sized, identifiable with standard thresholds
# Hard   = single-point mutations, small genes, or high sequence variability
CATEGORY_DIFFICULTY: dict[str, str] = {
    # VF categories
    "Exotoxin":           "Easy",    # large toxin genes (hla, lukF, sea, seb)
    "Motility":           "Easy",    # flagella = 7+ genes in one cluster
    "Secretion System":   "Easy",    # multi-gene T3SS / T6SS complexes
    "Adherence":          "Medium",  # adhesins vary in size
    "Biofilm":            "Medium",  # quorum sensing genes moderate size
    "Iron Uptake":        "Medium",  # siderophore biosynthesis clusters
    "Endotoxin":          "Medium",  # LPS genes, moderate
    "Invasion":           "Medium",  # invasins, moderate
    "Immune Evasion":     "Hard",    # capsule genes, high variability
    "Regulation":         "Hard",    # small regulator genes
    # AMR categories
    "Beta-Lactam Resistance":       "Easy",    # mecA, blaZ — large, well conserved
    "Efflux Pump":                  "Easy",    # large multi-drug efflux operons
    "Macrolide Resistance":         "Medium",  # erm genes, moderate
    "Tetracycline Resistance":      "Medium",  # tet genes, moderate
    "Glycopeptide Resistance":      "Medium",  # van cluster, moderate
    "Sulfonamide Resistance":       "Medium",  # sul genes, moderate
    "Aminoglycoside Resistance":    "Hard",    # aac/aph — often single mutation
    "Fluoroquinolone Resistance":   "Hard",    # point mutations in gyrA/parC
}

# ── TPM Normalization ────────────────────────────────────────────────────────
TPM_SCALING_FACTOR = 1_000_000

# ── Logging ───────────────────────────────────────────────────────────────────
LOG_LEVEL = os.getenv("NANOPORE_LOG_LEVEL", "INFO")
LOG_FORMAT = "%(asctime)s | %(name)-30s | %(levelname)-8s | %(message)s"
LOG_FILE = BASE_DIR / "pipeline.log"

# ── API ───────────────────────────────────────────────────────────────────────
API_HOST = "0.0.0.0"
API_PORT = 8000
