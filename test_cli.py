#!/usr/bin/env python3
"""Test CLI commands."""
import sys
sys.path.insert(0, '/mnt/temp2/kronckbm/gitrepos/VariantFeatures')

from click.testing import CliRunner
from variantfeatures.cli import main

runner = CliRunner()

# Test stats
print("=== Testing stats ===")
result = runner.invoke(main, ['stats'])
print(f"Exit code: {result.exit_code}")
print(result.output)
if result.exception:
    import traceback
    traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)

# Test query
print("\n=== Testing query -g KCNH2 ===")
result = runner.invoke(main, ['query', '-g', 'KCNH2'])
print(f"Exit code: {result.exit_code}")
print(result.output[:2000] if len(result.output) > 2000 else result.output)
if result.exception:
    import traceback
    traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)
