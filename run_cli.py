#!/usr/bin/env python3
"""Run CLI directly."""
import sys
import os

# Add parent to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from variantfeatures.cli import main
main()
