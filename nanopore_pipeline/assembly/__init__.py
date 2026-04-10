"""
Assembly module for Nanopore reads using Flye.
Provides de novo assembly to reduce overlapping read redundancy before alignment.
"""

from .flye_wrapper import FlyeAssembler, AssemblyResult

__all__ = ["FlyeAssembler", "AssemblyResult"]
