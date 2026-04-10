"""
Centralized logging configuration for the pipeline.
"""

import logging
import sys
from config.settings import LOG_LEVEL, LOG_FORMAT, LOG_FILE


def setup_logging(name: str = "nanopore_pipeline") -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    if not logger.handlers:
        # Console handler
        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(logging.Formatter(LOG_FORMAT))
        logger.addHandler(console)

        # File handler
        try:
            file_handler = logging.FileHandler(LOG_FILE, encoding="utf-8")
            file_handler.setFormatter(logging.Formatter(LOG_FORMAT))
            logger.addHandler(file_handler)
        except OSError:
            logger.warning("Could not create log file at %s", LOG_FILE)

    return logger
