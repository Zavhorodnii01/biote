"""
Database migration script to add assembly columns to existing Sample table.

Run this script to upgrade existing databases to support assembly metadata.
"""

import sqlite3
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from config import settings


def migrate_add_assembly_columns(db_path: Path):
    """Add assembly columns to Sample table if they don't exist."""

    print(f"Migrating database: {db_path}")

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Get existing columns
    cursor.execute("PRAGMA table_info(samples)")
    existing_columns = {row[1] for row in cursor.fetchall()}

    # Define new columns
    new_columns = [
        ("assembly_performed", "BOOLEAN DEFAULT 0"),
        ("assembly_path", "VARCHAR(1024)"),
        ("assembly_contigs", "INTEGER DEFAULT 0"),
        ("assembly_n50", "INTEGER"),
        ("assembly_total_bases", "INTEGER"),
        ("assembly_mean_contig_length", "FLOAT"),
        ("assembly_status", "VARCHAR(64)"),
    ]

    # Add missing columns
    added = 0
    for col_name, col_type in new_columns:
        if col_name not in existing_columns:
            try:
                sql = f"ALTER TABLE samples ADD COLUMN {col_name} {col_type}"
                cursor.execute(sql)
                print(f"  [+] Added column: {col_name}")
                added += 1
            except sqlite3.OperationalError as e:
                print(f"  [!] Error adding {col_name}: {e}")
        else:
            print(f"  [-] Column {col_name} already exists")

    conn.commit()
    conn.close()

    print(f"\nMigration complete: {added} columns added")


if __name__ == "__main__":
    db_path = Path(settings.DATABASE_URL.replace("sqlite:///", ""))

    if not db_path.exists():
        print(f"Database not found: {db_path}")
        print("Run this script after initializing the database.")
        sys.exit(1)

    migrate_add_assembly_columns(db_path)
