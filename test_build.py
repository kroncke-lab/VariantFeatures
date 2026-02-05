#!/usr/bin/env python3
"""Test the build command with gnomAD."""
import sys
sys.path.insert(0, '/mnt/temp2/kronckbm/gitrepos/VariantFeatures')

from click.testing import CliRunner
from variantfeatures.cli import main

runner = CliRunner()

print("=== Testing build -g KCNH2 --sources gnomad ===")
result = runner.invoke(main, ['build', '-g', 'KCNH2', '--sources', 'gnomad'])
print(f"Exit code: {result.exit_code}")
print(result.output)
if result.exception:
    import traceback
    traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)

# Check stats after build
print("\n=== Stats after build ===")
result = runner.invoke(main, ['stats'])
print(result.output)
