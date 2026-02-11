"""Catalytic triad/active site analysis module.

This module provides tools for analyzing catalytic triad and active site
geometry from MD trajectories.

Classes
-------
CatalyticTriadAnalyzer
    Main analyzer class for computing triad distances and contact fractions
"""

from polyzymd.analysis.triad.analyzer import CatalyticTriadAnalyzer

__all__ = ["CatalyticTriadAnalyzer"]
