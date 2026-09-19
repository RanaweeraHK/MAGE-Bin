"""Command-line interface for MAGE-Bin."""

from __future__ import annotations

import argparse
import importlib
import json
import shutil
import sys
from dataclasses import replace
from pathlib import Path

from . import __version__
from .config import MageBinConfig


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="magebin",
        description="Masked adaptive graph-evidence binning for viral metagenomes",
    )
    parser.add_argument(
        "--version", action="version", version=f"MAGE-Bin {__version__}"
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    bin_parser = subcommands.add_parser("bin", help="bin one processed cohort")
    bin_parser.add_argument("dataset", type=Path, help="model-ready dataset directory")
    bin_parser.add_argument("--output", "-o", type=Path, required=True)
    bin_parser.add_argument("--coverage", type=Path)

    run_parser = subcommands.add_parser(
        "run", help="prepare contigs, coverage and assembly graph, then bin"
    )
    run_parser.add_argument("--contigs", type=Path, required=True)
    run_parser.add_argument("--coverage", type=Path, required=True)
    run_parser.add_argument(
        "--graph", type=Path, required=True, help="GFA assembly graph"
    )
    run_parser.add_argument(
        "--paths", type=Path, help="SPAdes contigs.paths for segment-level GFA"
    )
    run_parser.add_argument("--output", "-o", type=Path, required=True)

    for command_parser in (bin_parser, run_parser):
        command_parser.add_argument(
            "--device", choices=("auto", "cpu", "cuda"), default="auto"
        )
        command_parser.add_argument("--seed", type=int, default=0)
        command_parser.add_argument("--min-contig-length", type=int, default=2_000)
        command_parser.add_argument("--minimum-epochs", type=int, default=40)
        command_parser.add_argument("--maximum-epochs", type=int, default=100)
        command_parser.add_argument("--graph-update-interval", type=int, default=10)
        command_parser.add_argument(
            "--no-complete-contig-gate",
            action="store_true",
            help="do not isolate contigs with terminal-repeat evidence",
        )

    doctor_parser = subcommands.add_parser(
        "doctor", help="report required Python packages and optional tools"
    )
    doctor_parser.add_argument("--json", action="store_true", dest="as_json")
    return parser


def _doctor(as_json: bool) -> int:
    required_modules = (
        "numpy",
        "pandas",
        "scipy",
        "sklearn",
        "torch",
        "igraph",
        "leidenalg",
    )
    optional_commands = ("checkv", "mmseqs", "samtools")
    packages: dict[str, dict[str, object]] = {}
    success = True
    for name in required_modules:
        try:
            module = importlib.import_module(name)
            packages[name] = {
                "available": True,
                "version": getattr(module, "__version__", "unknown"),
            }
        except (
            ImportError,
            OSError,
        ) as error:  # pragma: no cover - environment dependent
            packages[name] = {"available": False, "error": str(error)}
            success = False
    commands = {
        name: {"available": shutil.which(name) is not None, "path": shutil.which(name)}
        for name in optional_commands
    }
    report = {
        "magebin": __version__,
        "python": sys.version.split()[0],
        "required_python_packages": packages,
        "optional_commands": commands,
    }
    if as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(f"MAGE-Bin {__version__} | Python {report['python']}")
        for name, status in packages.items():
            detail = status.get("version", status.get("error", "missing"))
            marker = "ok" if status["available"] else "missing"
            print(f"  [{marker:7}] {name}: {detail}")
        print("Optional downstream commands:")
        for name, status in commands.items():
            marker = "ok" if status["available"] else "missing"
            print(f"  [{marker:7}] {name}: {status['path'] or '-'}")
    return 0 if success else 1


def main(argv: list[str] | None = None) -> int:
    """Run the MAGE-Bin command-line application."""

    parser = _build_parser()
    arguments = parser.parse_args(argv)
    if arguments.command == "doctor":
        return _doctor(arguments.as_json)
    if arguments.command in {"bin", "run"}:
        from .pipeline import run_binning

        config = replace(
            MageBinConfig(),
            seed=arguments.seed,
            min_contig_length=arguments.min_contig_length,
            minimum_epochs=arguments.minimum_epochs,
            maximum_epochs=arguments.maximum_epochs,
            graph_update_interval=arguments.graph_update_interval,
            isolate_complete_contigs=not arguments.no_complete_contig_gate,
        )
        if arguments.command == "run":
            from .preprocessing import prepare_dataset

            dataset = prepare_dataset(
                arguments.contigs,
                arguments.coverage,
                arguments.graph,
                arguments.output,
                min_contig_length=config.min_contig_length,
                paths=arguments.paths,
            )
            result = run_binning(
                dataset, arguments.output, config=config, device=arguments.device
            )
        else:
            result = run_binning(
                arguments.dataset,
                arguments.output,
                config=config,
                coverage_file=arguments.coverage,
                device=arguments.device,
            )
        print(
            f"MAGE-Bin assigned {result.contig_count} contigs to "
            f"{result.bin_count} bins ({result.singleton_count} singletons)."
        )
        print(f"Assignments: {result.assignment_file}")
        if result.bin_fasta_directory is not None:
            print(f"Bin FASTAs: {result.bin_fasta_directory}")
        print(f"Run metadata: {result.run_metadata_file}")
        return 0
    parser.error(f"unknown command: {arguments.command}")
    return 2
