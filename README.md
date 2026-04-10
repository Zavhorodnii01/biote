# Nanopore Pathogenicity Diagnostic Pipeline

**Thesis title (EN):** Detection of the pathogenicity in the clinical microbial Nanopore sequencing
**Thesis title (PL):** Wykrywanie patogenności w sekwencjonowaniu mikroorganizmów pochodzących z klinicznych próbek metodą Nanopore
**Status:** Active development — not final, still in discussion with supervisor.

---

## Meeting #1 — Supervisor Session

> This section is the working log for the thesis. Each sub-section maps to a timestamp from the recorded meeting. The next meeting is **Thursday, March 5** (Zoom, supervisor will send link).

---

### [0:09–1:29] What BLAST / Diamond Output Looks Like

**Summary:**
Alignment tools like BLAST and DIAMOND produce a tab-separated table. Each row is one hit: your query sequence matched against one reference entry. The file format is called "blast6" or "tab".

**Columns you get (already implemented in `AlignmentWrapper`):**

| Column | Meaning |
|--------|---------|
| `qseqid` | Your query ID (from Nanopore FASTA header) |
| `sseqid` | Reference ID (VFDB/CARD entry that matched) |
| `pident` | % identical positions in alignment |
| `length` | Alignment length |
| `mismatch` | Number of mismatches |
| `gapopen` | Number of gaps |
| `qstart/qend` | Where alignment starts/ends in your query |
| `sstart/send` | Where it starts/ends in the reference |
| `evalue` | Want as **small** as possible |
| `bitscore` | Want as **high** as possible |

**FASTA file format reminder:**
```
>sequence_header_1
ATCGATCGATCG...
>sequence_header_2
GCTAGCTAGCTA...
```

**Status:** Pipeline already handles this correctly via `AlignmentWrapper` using `blastn` and `diamond blastx`.

---

### [1:42–4:33] VFDB — Virulence Factor Database

**Summary:**
VFDB is a curated database of known virulence genes across pathogens. It is treated as frozen — download once, record the date, use that snapshot for the whole thesis.

**Rules from supervisor:**

- You **do not** need to keep updating it. Download it once.
- Write a **paragraph in your thesis** that includes:
  - What VFDB is and where it comes from (with proper citation)
  - Exact download date
  - How many sequences/proteins it contains
  - Which sections you kept and which you removed, and **why** (must be logical, e.g. remove protozoa if your work is only on bacteria)
- Size: small gigabytes — should be manageable
- Computational resources: supervisor will check availability at Politechnika, may run jobs for you if needed

**VFDB URLs (already in `config/settings.py`):**
- Nucleotide: `http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz`
- Protein: `http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz`

> **TODO (you):** Download VFDB. Record: exact download date, total sequences, list of sections present, which sections you removed and why. Write thesis paragraph.

---

### [4:49–5:41] Thesis Writing — Start Before Next Meeting

**Summary:**
No writing had been done yet. Supervisor expects a start before the March 5 meeting.

**What to do:**
- Set up **Overleaf** with thesis chapter/section structure (even skeleton headings are useful)
- Write **Introduction / Background** covering:
  - What Nanopore sequencing is and why it is introduced clinically
  - What pathogenicity means (virulence genes, AMR)
  - Why detecting virulence genes from sequencing matters for infectious disease diagnosis
- Supervisor will give feedback directly in Overleaf on what is missing

> **TODO (you):** Create Overleaf project. Write Introduction/Background draft. Bring to March 5 meeting.

---

### [5:50–12:58] Correct Pipeline Architecture & Alignment Method

**Summary:**
Student showed `pandas.merge()` being used on a simulated BLAST dataframe. Supervisor clarified that merge is only useful *after* you have real alignment output — it cannot replace running Diamond/BLAST. The correct pipeline flow was described:

**Correct pipeline order (what the thesis should implement):**

```
1. Input: Nanopore FASTA (raw long reads)
          ↓
2. Gene calling / prediction (find genes on the reads)
   → Tools: Prodigal, Pyrodigal, FragGeneScan
          ↓
3. Translate predicted genes → protein sequences
   (Diamond blastx can also handle translation internally)
          ↓
4. Align against VFDB protein sequences using Diamond blastx
   → Output: tabular file (qseqid, sseqid, pident, bitscore, etc.)
          ↓
5. Filter hits by identity %, bitscore, coverage
   → Keep only confident, real matches
          ↓
6. Classify surviving hits → virulence categories + risk score
          ↓
7. Report / visualise results
```

**On Diamond blastx vs blastn:**
- `blastx` = DNA query searched against protein DB. Diamond handles translation internally. **Recommended for VFDB** (protein sequences).
- `blastn` = DNA vs DNA. Used for CARD nucleotide models.

**Missing step (not yet in code):** Gene calling (step 2). Currently the pipeline sends raw reads directly to Diamond. This is an approximation — the more biologically correct approach is to predict genes first with Prodigal/Pyrodigal, then align the extracted proteins.

> **TODO (code, future):** Add gene calling step using Pyrodigal (Python bindings for Prodigal). Extract predicted protein sequences, then pipe those into `run_diamond_blastx()`.

---

### [12:58–17:00] Scientific Goal & LLM Component

**Summary:**
Supervisor outlined what the thesis scientific question could look like, and introduced a second, LLM-based pipeline track.

**Example scientific question:**
> "In gut microbiome samples, how has virulence content changed before vs. after COVID?"

This requires:
1. A pipeline that can quantify virulence genes from a Nanopore sample (what we are building)
2. Characterisation of virulence type (exotoxin, adhesion, immune evasion, etc.)
3. Cross-sample comparison (already partly implemented via `compare_samples()`)

**Two tracks being discussed:**

| Track | Approach | Status |
|-------|----------|--------|
| Traditional | Diamond → VFDB → classify | Implemented (gene calling step missing) |
| LLM-based | Feed sequence to LLM → predict virulence | To be explored |

**Website idea:**
A web UI where you paste a DNA or protein sequence and get back:
- Is this virulent? (Yes / No / Probability)
- What type of virulence?

Supervisor said this would be interesting specifically for the **LLM track**. The traditional pipeline can serve as a REST API backend for comparison.

> **TODO (you):** Research protein language models (ESM2, ProtTrans) and whether they can classify virulence from sequence. Bring ideas to March 5 meeting.

---

### [17:00–end] Alignment Quality Cutoffs — Justify Your Thresholds

**Summary:**
Any cutoff applied (minimum identity %, minimum bitscore) must be **justified with data**, not chosen arbitrarily. The way to do this is to plot score distributions and show where you drew the line.

**What supervisor described:**
> "You can plot, for example, the distribution and this line on the 90 and say, here, see, I'm cutting off this tail and this is real stuff."

**Concept: best-hit selection**
A single query gene can match multiple VFDB entries. You need to pick the best one:
- 100% identity + full-length coverage → definitely the same gene
- 60% identity, 80% coverage → uncertain, may or may not be real
- Use bitscore to rank and keep only the top match per query

**For the thesis:**
1. Run alignment → get all raw hits
2. Plot histogram of `pident` (% identity) across all hits
3. Plot histogram of `bitscore` across all hits
4. Choose cutoffs based on what you see (natural break in distribution)
5. Draw the cutoff line on the histogram
6. State in thesis: "I removed all hits below X% identity because the distribution shows a natural boundary at this point"

---

## Task List

### Your Tasks (manual — cannot be automated)

- [ ] Download VFDB — record download date, sequence count, sections
- [ ] Write VFDB thesis paragraph (citation, date, sections selected/removed + reason)
- [ ] Set up Overleaf with thesis structure
- [ ] Write Introduction / Background section (bring to March 5)
- [ ] Research LLM approaches for virulence classification from sequence
- [ ] Read about alignment statistics in more depth (pident, bitscore, coverage, evalue)

### Done Automatically (code changes applied while you slept)

- [x] **`select_best_hits()`** added to `AlignmentWrapper`
  For each query gene, only the best VFDB hit (by bitscore) is kept. Addresses supervisor's question: "you will have the same gene aligned against seven different hits — how would you know which one is good?"
  → `nanopore_pipeline/alignment/wrapper.py`

- [x] **`plot_alignment_quality()`** added to `PathogenicityReporter`
  Generates two interactive histograms per sample run:
  - Distribution of `% identity` across all hits (with adjustable cutoff line)
  - Distribution of `bitscore` across all hits (with adjustable cutoff line)
  These charts are saved to `data/results/` and can be used directly as thesis figures.
  Output: `{sample}_vf_alignment_quality.html` and `{sample}_amr_alignment_quality.html`
  → `nanopore_pipeline/reporter/visualisation.py`

- [x] **Pipeline runner updated** to generate quality charts automatically after each alignment step
  → `nanopore_pipeline/pipeline_runner.py`

- [x] **API endpoint updated** to handle new return type from `align_sample()`
  → `nanopore_pipeline/api/endpoints.py`

---

## Research: Methods for Virulence Classification from Sequences

> Literature review: *"Evolutionary and Computational Architectures for the Classification of Microbial Virulence from Biological Sequences"*
> This table summarises every approach covered in the research and marks what is directly useful for your thesis.

---

### Track 1 — Homology-Based (Traditional, Already Implemented)

| Tool | Input | What it does | Status in this project |
|------|-------|-------------|----------------------|
| BLAST (`blastn`) | DNA vs DNA | Aligns Nanopore reads against nucleotide DB (CARD) | **Done** |
| Diamond (`blastx`) | DNA vs Protein | Aligns reads against VFDB protein sequences, handles translation | **Done** |

**How it works:** Looks for sequences similar to known virulence genes. Fails on novel pathogens with no close relatives in the database.

---

### Track 2 — Genomic Language Models (DNA-level)

| Model | Tokenization | Key strength | Weakness | Relevance |
|-------|-------------|-------------|----------|-----------|
| **DNABERT** | Overlapping k-mer (k=6) | Good at local motifs (promoters, splice sites) | High memory, slow (O(n²)) | Low — too slow for long Nanopore reads |
| **DNABERT-2** | Byte Pair Encoding (BPE) | 21x fewer params than SOTA, multi-species, handles long seqs | Still limited context | **Medium — good for DNA-level virulence detection** |
| **HyenaDNA** | Nucleotide-level | Up to 1 million token context, sub-quadratic (O(n log n)) | Less interpretable | **Medium — best for long Nanopore reads** |
| **Nucleotide Transformer** | Non-overlapping k-mer | Cross-species generalisation, massive pre-training | Large model size | Low — overkill for this scope |

**What BPE means for you:** Instead of splitting `ATCGATCG` into overlapping `ATCGAT`, `TCGATC`, `CGATCG`... it merges common patterns into single tokens. Result: ~4-5x shorter input → faster processing.

---

### Track 3 — Protein Language Models (Protein-level)

| Model | Architecture | Key strength | Best for |
|-------|-------------|-------------|---------|
| **ESM-2 (650M params)** | 33-layer Transformer | Learns structural + evolutionary context per amino acid | Gram-positive bacteria virulence |
| **ESM-1b** | Transformer | Strong residue-residue interaction capture | Gram-negative bacteria virulence |
| **ProtBert / ProtBert-BFD** | BERT masked LM | Good general baseline, widely used | Either, good starting point |
| **ESMFold** | ESM-2 + folding head | Predicts 3D structure from sequence | Structure-aware virulence detection |

**Key finding from research:** G+ and G- bacteria need different models.
- G+ (e.g. *Staphylococcus*, *Streptococcus*) → ESM-2-650M
- G- (e.g. *E. coli*, *Klebsiella*) → ESM-1b

---

### Track 3a — Ready-to-Use Virulence-Specific Tools

These are models already fine-tuned specifically for virulence prediction — you don't train them from scratch.

| Tool | Input | Approach | Accuracy | Relevant? |
|------|-------|----------|----------|-----------|
| **pLM4VF** | Protein sequence | ESM-2/ESM-1b fine-tuned, separate G+/G- models | 76.2% (G+), 82.2% (G-) | **High — plug-and-play for your LLM track** |
| **PLMVF** | Protein sequence + 3D structure | ESM-2 embeddings + ESMFold structure + KAN classifier | 86.1% | **High — most accurate, handles remote homology** |
| **DTVF** | Protein sequence | ProtT5 + CNN/LSTM transfer learning | 84.55% | Medium — good accuracy, easier than PLMVF |
| **VirulentPred 2.0** | Protein sequence + PSSM | Ensemble ML with position-specific scoring | 85.18% | Medium — older approach, PSSM requires alignment |
| **SeqScreen** | Protein sequence | Ensemble ML → FunSoC labels (cytotoxicity, secretion, etc.) | — | **High — gives functional category labels, not just binary** |

**Remote homology (PLMVF specialty):** Can detect virulence in proteins with <25% sequence identity to known virulence factors — this is exactly what alignment-based BLAST misses.

---

### Track 4 — AMR-Specific Models

| Tool | Input | Approach | Notes |
|------|-------|----------|-------|
| **ProtBert-BFD + MH-LSTM** | Full proteome | Per-protein embeddings summed for phenotypic resistance prediction | Good for strain-level AMR |
| **Evo-MoE** | Genome sequence | Mixture of Experts + genetic algorithm simulates mutation/selection | Predicts evolutionary trajectories of resistance — research-level |

---

### Comparison: Traditional vs LLM Track

| Aspect | BLAST/Diamond (Track 1) | LLM-based (Track 2/3) |
|--------|------------------------|-----------------------|
| Detects known virulence genes | Yes | Yes |
| Detects novel/divergent pathogens | No | Yes |
| Requires reference database | Yes (VFDB) | No (pre-trained on UniRef) |
| Speed | Fast | Slow (GPU needed) |
| Interpretability | High (you see the match) | Lower (embedding space) |
| Accuracy on known genes | High | Comparable |
| Accuracy on novel genes | Fails | Works |
| Thesis novelty | Low | **High** |

---

### What to Actually Use in Your Thesis

| Priority | Tool | Why | How hard |
|----------|------|-----|----------|
| **Done** | BLAST + Diamond + VFDB | Track 1 baseline, already works | Done |
| **Do next** | pLM4VF | Fine-tuned virulence predictor, plug-and-play, compare vs BLAST | Medium — pip install + inference |
| **Do next** | SeqScreen | Gives FunSoC labels (cytotoxic, secretion system, etc.) — good for your category system | Medium |
| **Stretch goal** | PLMVF | Best accuracy (86.1%), handles novel pathogens, good thesis story | Hard — needs ESMFold |
| **For the website** | ESM-2 + fine-tune or pLM4VF | User pastes sequence → virulence probability | Medium |
| **Skip for now** | DNABERT, HyenaDNA | DNA-level, need gene calling first | Complex |
| **Skip for now** | Evo-MoE | AMR evolution simulation — too far from core scope | Very hard |

---

### Key Concepts to Know for Your Thesis

| Concept | What it means |
|---------|--------------|
| **Tokenization** | How sequences are split into pieces the model processes |
| **k-mer** | Fixed-length overlapping substrings (e.g. all 6-mers of a sequence) |
| **BPE** | Smarter tokenization — merges frequent patterns, shorter input = faster |
| **Embedding** | A vector (list of numbers) representing a sequence in the model's learned space |
| **Remote homology** | Two proteins with the same function but <25% sequence similarity — BLAST misses these, pLMs find them |
| **FunSoC** | Functional label: "cytotoxicity", "secretion", "immune evasion" — more useful than just "virulent/not" |
| **Fine-tuning** | Taking a pre-trained model and training it a bit more on virulence-specific data |
| **Zero-shot** | Using the model as-is without any fine-tuning — possible with large enough models |
| **TPM** | Normalisation for cross-sample comparison (already in your pipeline) |

---

## Next Meeting

**Date:** Thursday, March 5
**Format:** Zoom (supervisor will send link)
**Bring:**
- Overleaf link with thesis structure and introduction draft
- VFDB download info (date, sections, size)
- LLM ideas for the sequence classification track

---

## Key Rules & Notes

| Topic | Rule |
|-------|------|
| VFDB | Frozen at download date. Cite properly. Record date and sequence count. |
| Cutoffs | Never arbitrary — plot distributions, justify with data, show cutoff line on histogram. |
| Best-hit | Per query gene, keep only the single best VFDB match by bitscore. |
| Diamond blastx | Preferred for VFDB — DNA query vs protein DB, Diamond handles translation. |
| `--long-reads` | Required for Nanopore reads in Diamond (already set). |
| Gene calling | Should come before alignment (Prodigal/Pyrodigal). Not yet implemented. |
| E-value | Lower is better. Bitscore: higher is better. |
| Writing | Start now. Even skeleton Overleaf structure is useful for supervisor feedback. |
| LLM track | More novel than traditional track — interesting for a website/comparison study. |
| TPM | Already implemented — enables fair cross-sample comparison. |

---

## Pipeline — Current Outputs Per Run

```
data/results/
├── {sample}_profile.html                   # Hit count + TPM bar chart by category
├── {sample}_risk_dashboard.html            # Pie chart (risk distribution) + detail table
├── {sample}_vf_alignment_quality.html      # NEW — pident & bitscore histograms (VF)
├── {sample}_amr_alignment_quality.html     # NEW — pident & bitscore histograms (AMR)
└── {sample}_classification.json            # JSON export of all category profiles
```

---

A production-grade, full-stack bioinformatics pipeline for classifying Oxford Nanopore sequencing results as clinically relevant pathogens. The tool doesn't just identify a species -- it evaluates pathogenic potential by analysing functional gene profiles against **VFDB** (Virulence Factor Database) and **CARD** (Comprehensive Antibiotic Resistance Database).

**Input format: FASTA** (basecalled Nanopore reads).

---

## Initial Prompt

> Develop a pipeline that takes Nanopore long-read sequencing data (FASTA), aligns reads against VFDB and CARD reference databases, classifies hits into pathogenicity categories (Adherence, Motility, Exotoxins, Endotoxins, Secretion Systems, AMR), computes TPM-normalised abundance, and provides comparative analysis with interactive visualisations via a REST API.

## What Was Built

### Full Pipeline Flow

```
Raw FASTA ──> Read Stats ──> Register Sample ──> BLAST vs VFDB ──> BLAST vs CARD ──> Classify ──> Report
   │              │                │                   │                  │              │            │
   │         total_reads      store in DB         VF hits to DB     AMR hits to DB   TPM + Risk   HTML + JSON
   │         N50, GC%                                                                  levels      charts
   │         mean_length
   └── Single endpoint: POST /run (FASTA in -> full result out)
```

### Architecture

```
Bioinformatics/
├── main.py                          # FastAPI entry point (uvicorn)
├── requirements.txt                 # Python dependencies
├── config/
│   └── settings.py                  # Central configuration (paths, thresholds, categories)
├── nanopore_pipeline/
│   ├── pipeline_runner.py           # End-to-end orchestrator (FASTA in -> report out)
│   ├── models/
│   │   └── database.py              # SQLAlchemy ORM (5 tables: VFs, AMR, Samples, Hits, Classifications)
│   ├── db/
│   │   └── manager.py               # Module 1: Database Manager (download/parse/index VFDB & CARD)
│   ├── alignment/
│   │   └── wrapper.py               # Module 2: BLAST+/DIAMOND wrapper (Nanopore-optimised params)
│   ├── classifier/
│   │   └── pathogenicity.py         # Module 3: Pathogenicity Classifier (TPM, risk scoring, comparison)
│   ├── reporter/
│   │   └── visualisation.py         # Module 4: Plotly-based visualisation suite
│   ├── api/
│   │   └── endpoints.py             # FastAPI REST endpoints (including /run for full pipeline)
│   ├── utils/
│   │   ├── logging_config.py        # Centralised logging
│   │   └── fasta_parser.py          # FASTA parser with length filtering, N50/GC stats
│   └── tests/
│       ├── test_fasta_parser.py     # Parser + stats unit tests
│       ├── test_database.py         # DB model & manager tests
│       ├── test_classifier.py       # Classifier + TPM tests
│       ├── test_alignment.py        # BLAST output parsing & filtering tests
│       └── test_pipeline_runner.py  # End-to-end pipeline tests
└── data/
    ├── vfdb/                        # VFDB reference files
    ├── card/                        # CARD reference files
    ├── samples/                     # Uploaded Nanopore FASTA files
    └── results/                     # Generated reports (HTML, JSON)
```

### Module Details

#### Pipeline Runner (`nanopore_pipeline/pipeline_runner.py`)
- Single entry point: `PipelineRunner.run(sample_name, fasta_path)` -> `PipelineResult`
- Chains all 6 steps: stats -> register -> BLAST VFDB -> BLAST CARD -> classify -> report
- Also exposes `run_comparison()` for two-sample analysis
- Can skip VF or AMR alignment independently via flags

#### Module 1: Database Manager (`nanopore_pipeline/db/manager.py`)
- Downloads and decompresses VFDB/CARD reference FASTA files
- Parses FASTA headers to extract gene names, organisms, and functional annotations
- Auto-categorises virulence factors into 10 categories (Adherence, Motility, Exotoxin, Endotoxin, Secretion Systems, Iron Uptake, Immune Evasion, Biofilm, Invasion, Regulation)
- Auto-categorises AMR genes into 8 drug-class categories
- Idempotent loading (no duplicates on re-run)
- SQLite by default, PostgreSQL-ready via `DATABASE_URL` env var

#### Module 2: Alignment Wrapper (`nanopore_pipeline/alignment/wrapper.py`)
- Python interface for BLAST+ (`blastn`) and DIAMOND (`blastx`)
- **Nanopore-optimised parameters**: word_size=7, reward=1, penalty=-1, gapopen=2, gapextend=1 (tuned for ~10% error rate)
- BLAST database creation (`makeblastdb` / `diamond makedb`)
- Tabular output parsing (outfmt 6 with qlen/slen)
- Configurable filtering by percent identity and query coverage
- Results stored as ORM objects linked to samples and reference entries

#### Module 3: Pathogenicity Classifier (`nanopore_pipeline/classifier/pathogenicity.py`)
- Groups alignment hits by VF category and AMR drug class
- **TPM (Transcripts Per Million) normalisation** for cross-sample comparability
- **Risk scoring**:
  - `Critical` -- Multiple critical categories present (Exotoxin + Secretion System + Beta-Lactam)
  - `High` -- One critical or multiple high-risk categories
  - `Medium` -- Moderate hit counts in non-critical categories
  - `Low` -- Minimal or non-specific hits
- **Comparative analysis**: side-by-side TPM fold-change between two samples
- Results persisted to `classification_results` table

#### Module 4: Comparative Reporter (`nanopore_pipeline/reporter/visualisation.py`)
- Interactive Plotly charts:
  - Category bar charts (hit counts + TPM side by side)
  - Comparative heatmaps (sample A vs. sample B)
  - Fold-change bar charts (log-scale, colour-coded)
  - Risk dashboard (pie chart + detail table)
- JSON export for downstream integration
- All charts saved as standalone HTML files

#### FASTA Parser (`nanopore_pipeline/utils/fasta_parser.py`)
- Parses plain and gzipped FASTA files
- Computes per-file statistics: total reads, total bases, mean read length, N50, GC content
- Length-based read filtering (default min 200 bp for Nanopore)

### REST API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| `POST` | `/run` | **Full pipeline**: upload FASTA -> stats + align + classify + report |
| `GET` | `/` | Health check |
| `POST` | `/db/init` | Initialise database |
| `POST` | `/db/load-vfdb` | Load VFDB reference data |
| `POST` | `/db/load-card` | Load CARD reference data |
| `POST` | `/samples/upload` | Upload a FASTA file (register only) |
| `GET` | `/samples` | List all samples |
| `GET` | `/samples/{name}` | Get sample details |
| `POST` | `/align` | Run BLAST/DIAMOND alignment (step-by-step) |
| `POST` | `/classify/{name}` | Classify a sample (step-by-step) |
| `GET` | `/compare?sample_a=X&sample_b=Y` | Compare two samples |
| `GET` | `/report/{name}` | Generate interactive HTML report |
| `GET` | `/report/compare/{a}/{b}` | Generate comparison HTML report |

### Nanopore-Specific Design Decisions
- **Input format**: FASTA (basecalled Nanopore output, no quality scores)
- **Read filtering**: minimum length (200 bp) based on sequence length only
- **Read statistics**: N50, GC content, mean/min/max read length computed from FASTA
- **BLAST parameters**: small word size and permissive gap penalties for high-error long reads
- **DIAMOND `--long-reads` flag**: enables Nanopore-aware alignment mode
- **No dust/low-complexity filtering**: disabled to avoid masking real Nanopore signal

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run the API server
python main.py
# or: uvicorn main:app --reload

# Open Swagger UI
# http://localhost:8000/docs

# Run tests
pytest nanopore_pipeline/tests/ -v
```

### Programmatic Usage

```python
from pathlib import Path
from nanopore_pipeline.pipeline_runner import PipelineRunner

runner = PipelineRunner()

# Run full pipeline on a FASTA file
result = runner.run("patient_01", Path("reads.nucleotide_fasta_protein_homolog_model.fasta"))
print(f"Reads: {result.fasta_stats.total_reads}")
print(f"N50: {result.fasta_stats.n50}")
print(f"VF hits: {result.vf_hits}, AMR hits: {result.amr_hits}")
for p in result.profiles:
    print(f"  {p.category}: {p.hit_count} hits, risk={p.risk_level}")

# Compare two samples
comparison = runner.run_comparison("patient_01", "patient_02")
```

## External Tool Requirements
- **BLAST+** (ncbi-blast) -- for `blastn` / `makeblastdb`
- **DIAMOND** -- for protein-level alignment (optional)
- Both are expected to be on `$PATH`. The pipeline works without them for classification/reporting on pre-computed hits.

## Next Steps (for future prompts)
- Assembly integration (Flye for long-read assembly)
- Prokka annotation layer
- Streamlit dashboard as alternative UI
- Docker containerisation
- Batch processing for multiple samples
- PDF report generation
