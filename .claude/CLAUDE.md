# Bioinformatics Thesis Project - Context

## Student Profile
- **Role**: Medical Informatics Student (Politechnika)
- **Project**: Diploma Thesis - Nanopore Pathogenicity Detection
- **Thesis Title (EN)**: Detection of the pathogenicity in the clinical microbial Nanopore sequencing
- **Thesis Title (PL)**: Wykrywanie patogenności w sekwencjonowaniu mikroorganizmów pochodzących z klinicznych próbek metodą Nanopore
- **Status**: Active development - Work in progress, not final
- **Next Supervisor Meeting**: Check README.md for latest date

## Project Overview
This is a **production-grade bioinformatics pipeline** that analyzes Oxford Nanopore long-read sequencing data to detect pathogenic potential in clinical microbial samples. The pipeline doesn't just identify species - it evaluates pathogenic potential by analyzing functional gene profiles.

### What the Pipeline Does
```
FASTA reads → Stats → Align vs VFDB → Align vs CARD → Classify → Report (HTML + JSON)
     ↓           ↓           ↓              ↓            ↓            ↓
  Raw data   N50, GC%   Virulence    Resistance    TPM norm'd   Interactive
                        genes        genes         risk scores   visualizations
```

### Input/Output
- **Input**: FASTA files (basecalled Nanopore reads)
- **Output**:
  - Interactive HTML reports (Plotly charts)
  - JSON exports
  - Risk classifications (Critical/High/Medium/Low)
  - TPM-normalized gene abundance
  - Quality histograms (% identity, bitscore)

## Key Technologies & Databases

### Reference Databases
1. **VFDB (Virulence Factor Database)** - Curated virulence genes
   - Frozen at download date (record in thesis)
   - Categories: Adherence, Motility, Exotoxin, Endotoxin, Secretion Systems, etc.

2. **CARD (Comprehensive Antibiotic Resistance Database)** - AMR genes
   - Drug class categorization (Beta-Lactam, Fluoroquinolone, etc.)

### Alignment Tools
- **BLAST+ (blastn)** - DNA vs DNA alignment
- **DIAMOND (blastx)** - DNA vs Protein alignment (with translation)
- Nanopore-optimized parameters: word_size=7, --long-reads flag

### Core Concepts
- **TPM (Transcripts Per Million)** - Normalization for cross-sample comparison
- **N50** - Read length quality metric
- **Risk Scoring** - Multi-factor pathogenicity assessment
- **Best-hit selection** - Per query gene, keep only top match by bitscore

## Architecture

```
Bioinformatics/
├── main.py                          # FastAPI server entry point
├── run_comparison.py                # 3-way comparison script (perfect/simulated/real)
├── simulate_reads.py                # Nanopore read simulator
├── prepare_samples.py               # Data preparation utilities
├── config/
│   └── settings.py                  # Central configuration
├── nanopore_pipeline/
│   ├── pipeline_runner.py           # Orchestrator (6-step flow)
│   ├── models/database.py           # SQLAlchemy ORM (5 tables)
│   ├── db/manager.py                # Database manager (download/parse/index)
│   ├── alignment/wrapper.py         # BLAST/DIAMOND wrapper
│   ├── classifier/pathogenicity.py  # TPM normalization + risk scoring
│   ├── reporter/visualisation.py    # Plotly charts + JSON export
│   ├── api/endpoints.py             # REST API endpoints
│   ├── utils/
│   │   ├── fasta_parser.py          # FASTA stats (N50, GC%, etc.)
│   │   └── logging_config.py
│   └── tests/                       # Pytest test suite
└── data/
    ├── vfdb/                        # VFDB reference files
    ├── card/                        # CARD reference files
    ├── samples/                     # Input FASTA files
    └── results/                     # Generated reports
```

## Pipeline Flow (6 Steps)

1. **Compute FASTA Stats** - total_reads, N50, GC%, mean_length
2. **Register Sample** - Store metadata in SQLite DB
3. **Align vs VFDB** - Find virulence factor matches (blastn/diamond)
4. **Align vs CARD** - Find AMR gene matches
5. **Classify** - TPM normalization + risk level assignment
6. **Report** - Generate interactive HTML charts + JSON

## Comparison Script (run_comparison.py)

**Purpose**: 3-way performance test comparing:
1. **Test 1**: Perfect reference DNA (no errors)
2. **Test 2**: Simulated Nanopore reads (~10% error)
3. **Test 3**: Real SRA Nanopore data (E. coli O157:H7)

**What it does**:
- Runs full pipeline on all three inputs
- Generates side-by-side comparison table
- Shows how alignment quality affects hit detection
- Uses different BLAST word_size per test (11 for perfect, 7 for Nanopore)

**Output**: Console table + HTML reports for each sample

## Research Tracks

### Track 1: Traditional Homology (IMPLEMENTED)
- BLAST/DIAMOND alignment against VFDB/CARD
- Fast, interpretable, works well for known pathogens
- **Missing**: Gene calling step (Prodigal/Pyrodigal) - should predict genes first

### Track 2: LLM-Based (FUTURE WORK)
Protein language models for virulence classification:
- **pLM4VF** - Fine-tuned ESM-2, separate G+/G- models
- **PLMVF** - ESM-2 + ESMFold + KAN classifier (86% accuracy)
- **SeqScreen** - Gives functional labels (cytotoxic, secretion, etc.)
- **DNABERT-2** - DNA-level LLM with BPE tokenization

**Advantage**: Can detect novel pathogens with low sequence similarity (remote homology)

### Track 3: Web Interface (PLANNED)
User pastes sequence → virulence probability + functional category

## Important Thesis Rules (from Supervisor)

1. **VFDB must be frozen** - Download once, record date, justify which sections you kept/removed
2. **Cutoffs must be justified with data** - Plot distributions, show natural breaks
3. **Best-hit selection** - Per query gene, keep only the top VFDB match
4. **Start writing NOW** - Set up Overleaf, write Introduction/Background
5. **Gene calling missing** - Should use Prodigal before alignment (not yet implemented)
6. **Alignment quality plots** - Use generated histograms (pident, bitscore) in thesis

## TODOs for Student (Manual)

- [ ] Download VFDB - record date, sequence count, sections
- [ ] Write VFDB thesis paragraph (citation, justification)
- [ ] Set up Overleaf with thesis structure
- [ ] Write Introduction/Background section
- [ ] Research LLM approaches (pLM4VF, PLMVF, SeqScreen)
- [ ] Prepare for next supervisor meeting

## API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/run` | **Full pipeline**: FASTA → complete result |
| POST | `/db/init` | Initialize database |
| POST | `/db/load-vfdb` | Load VFDB reference |
| POST | `/db/load-card` | Load CARD reference |
| POST | `/samples/upload` | Upload FASTA (register only) |
| GET | `/samples` | List all samples |
| POST | `/align` | Run alignment (step-by-step) |
| POST | `/classify/{name}` | Classify sample |
| GET | `/compare?sample_a=X&sample_b=Y` | Compare two samples |
| GET | `/report/{name}` | Generate HTML report |

## Key Files to Know

- **README.md** - Comprehensive thesis documentation (meeting notes, research, TODOs)
- **main.py** - FastAPI server entry
- **run_comparison.py** - 3-way test comparison
- **pipeline_runner.py** - Main orchestrator (6-step flow)
- **config/settings.py** - All configuration (paths, thresholds, categories)

## Common Operations

### Run the API server
```bash
python main.py
# or: uvicorn main:app --reload
# Access: http://localhost:8000/docs (Swagger UI)
```

### Run 3-way comparison
```bash
python run_comparison.py
```

### Run tests
```bash
pytest nanopore_pipeline/tests/ -v
```

### Programmatic usage
```python
from nanopore_pipeline.pipeline_runner import PipelineRunner
runner = PipelineRunner()
result = runner.run("sample_name", Path("reads.fasta"))
```

## Future Skills to Develop

### 1. Bioinformatics Database Fetcher
- Fetch sequences from NCBI, UniProt, VFDB, CARD
- Download and parse reference genomes
- Handle API rate limits and authentication

### 2. LLM Virulence Predictor
- Integrate pLM4VF or PLMVF models
- Sequence → virulence probability + functional category
- Compare LLM predictions vs BLAST results

### 3. Gene Caller Integration
- Wrap Prodigal/Pyrodigal for gene prediction
- Extract predicted proteins before alignment
- More biologically correct than aligning raw reads

### 4. Report Generator
- Generate LaTeX thesis sections automatically
- Create publication-ready figures from Plotly charts
- Statistical summaries for Methods section

## Communication Style
- This is thesis work - expect iterative development
- Not everything is final - plans may change after supervisor meetings
- Focus on justifying decisions with data (plots, statistics)
- Prioritize reproducibility and documentation

## Context Awareness
- The student may forget what they implemented previously
- Always check README.md for latest supervisor instructions
- Reference specific file locations when discussing code
- Explain bioinformatics concepts when relevant (TPM, N50, remote homology, etc.)
