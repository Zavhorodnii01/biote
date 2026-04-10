"""
Verification script for assembly integration.
Tests that all components are properly connected without requiring Flye.
"""

import sys
from pathlib import Path

print("=" * 70)
print("Assembly Integration Verification")
print("=" * 70)
print()

# Test 1: Import all modules
print("Test 1: Module imports...")
try:
    from nanopore_pipeline.assembly import FlyeAssembler, AssemblyResult
    from nanopore_pipeline.pipeline_runner import PipelineRunner
    from nanopore_pipeline.db.manager import DatabaseManager
    from config import settings
    print("  [OK] All modules imported successfully")
except ImportError as e:
    print(f"  [FAIL] Import failed: {e}")
    sys.exit(1)

# Test 2: Check configuration
print("\nTest 2: Configuration settings...")
try:
    assert hasattr(settings, 'ASSEMBLIES_DIR')
    assert hasattr(settings, 'FLYE_GENOME_SIZE')
    assert hasattr(settings, 'FLYE_MIN_OVERLAP')
    assert hasattr(settings, 'FLYE_ITERATIONS')
    assert hasattr(settings, 'FLYE_THREADS')
    assert hasattr(settings, 'FLYE_META_MODE')
    assert hasattr(settings, 'FLYE_MIN_READ_LENGTH')
    assert hasattr(settings, 'FLYE_TIMEOUT')
    print(f"  [OK] ASSEMBLIES_DIR: {settings.ASSEMBLIES_DIR}")
    print(f"  [OK] FLYE_GENOME_SIZE: {settings.FLYE_GENOME_SIZE}")
    print(f"  [OK] FLYE_MIN_OVERLAP: {settings.FLYE_MIN_OVERLAP}")
    print(f"  [OK] FLYE_TIMEOUT: {settings.FLYE_TIMEOUT}s")
except (AttributeError, AssertionError) as e:
    print(f"  [FAIL] Configuration check failed: {e}")
    sys.exit(1)

# Test 3: Check assemblies directory
print("\nTest 3: Assemblies directory...")
if settings.ASSEMBLIES_DIR.exists():
    print(f"  [OK] Directory exists: {settings.ASSEMBLIES_DIR}")
else:
    print(f"  [INFO] Creating directory: {settings.ASSEMBLIES_DIR}")
    settings.ASSEMBLIES_DIR.mkdir(parents=True, exist_ok=True)
    print(f"  [OK] Directory created")

# Test 4: Database schema
print("\nTest 4: Database schema...")
try:
    from nanopore_pipeline.models.database import Sample, init_db

    # Check Sample model has new fields
    sample_columns = Sample.__table__.columns.keys()
    required_columns = [
        'assembly_performed',
        'assembly_path',
        'assembly_contigs',
        'assembly_n50',
        'assembly_total_bases',
        'assembly_mean_contig_length',
        'assembly_status',
    ]

    for col in required_columns:
        if col in sample_columns:
            print(f"  [OK] Column '{col}' exists in Sample table")
        else:
            print(f"  [FAIL] Column '{col}' MISSING from Sample table")
            print(f"    Run: python nanopore_pipeline/db/migrate_add_assembly_columns.py")
            sys.exit(1)

except Exception as e:
    print(f"  [FAIL] Database schema check failed: {e}")
    sys.exit(1)

# Test 5: FlyeAssembler initialization
print("\nTest 5: FlyeAssembler initialization...")
try:
    assembler = FlyeAssembler()
    print(f"  [OK] FlyeAssembler created")
    print(f"    genome_size={assembler.genome_size}")
    print(f"    min_overlap={assembler.min_overlap}")
    print(f"    meta_mode={assembler.meta_mode}")
except Exception as e:
    print(f"  [FAIL] FlyeAssembler initialization failed: {e}")
    sys.exit(1)

# Test 6: PipelineRunner with assembly parameter
print("\nTest 6: PipelineRunner integration...")
try:
    runner = PipelineRunner()

    # Check that run() method accepts perform_assembly parameter
    import inspect
    sig = inspect.signature(runner.run)
    params = sig.parameters

    if 'perform_assembly' in params:
        print(f"  [OK] PipelineRunner.run() has 'perform_assembly' parameter")
        default_value = params['perform_assembly'].default
        print(f"    Default value: {default_value}")
    else:
        print(f"  [FAIL] PipelineRunner.run() missing 'perform_assembly' parameter")
        sys.exit(1)

except Exception as e:
    print(f"  [FAIL] PipelineRunner check failed: {e}")
    sys.exit(1)

# Test 7: DatabaseManager register_sample signature
print("\nTest 7: DatabaseManager.register_sample() signature...")
try:
    db = DatabaseManager()

    import inspect
    sig = inspect.signature(db.register_sample)
    params = sig.parameters

    required_params = [
        'assembly_performed',
        'assembly_path',
        'assembly_contigs',
        'assembly_n50',
        'assembly_total_bases',
        'assembly_mean_contig_length',
        'assembly_status',
    ]

    for param in required_params:
        if param in params:
            print(f"  [OK] Parameter '{param}' present")
        else:
            print(f"  [FAIL] Parameter '{param}' MISSING")
            sys.exit(1)

except Exception as e:
    print(f"  [FAIL] DatabaseManager check failed: {e}")
    sys.exit(1)

# Test 8: API endpoint check
print("\nTest 8: API endpoint integration...")
try:
    from nanopore_pipeline.api.endpoints import app

    # Get the /run endpoint
    run_route = None
    for route in app.routes:
        if hasattr(route, 'path') and route.path == '/run':
            run_route = route
            break

    if run_route:
        print(f"  [OK] /run endpoint exists")

        # Check endpoint signature
        import inspect
        endpoint_func = run_route.endpoint
        sig = inspect.signature(endpoint_func)

        if 'perform_assembly' in sig.parameters:
            print(f"  [OK] /run endpoint has 'perform_assembly' parameter")
        else:
            print(f"  [FAIL] /run endpoint missing 'perform_assembly' parameter")
    else:
        print(f"  [FAIL] /run endpoint not found")

except Exception as e:
    print(f"  [INFO] API check skipped (optional): {e}")

# Test 9: Check Flye availability (optional)
print("\nTest 9: Flye availability check...")
try:
    import shutil
    if shutil.which("flye"):
        print(f"  [OK] Flye is installed and available")
    else:
        print(f"  [INFO] Flye not found in PATH")
        print(f"    Assembly will fail at runtime without Flye")
        print(f"    Install with: conda install -c bioconda flye")
except Exception as e:
    print(f"  [INFO] Flye check failed: {e}")

# Final summary
print("\n" + "=" * 70)
print("VERIFICATION COMPLETE")
print("=" * 70)
print()
print("[OK] All critical components are properly integrated")
print()
print("Next steps:")
print("  1. Install Flye: conda install -c bioconda flye")
print("  2. Run tests: pytest nanopore_pipeline/tests/test_assembly.py -v")
print("  3. Try demo: python demo_assembly.py --help")
print()
