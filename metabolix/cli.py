"""Command line interface for the Metabolix BIDS application."""

from __future__ import annotations

import argparse
import logging
import os
import sys
import textwrap
from pathlib import Path
from typing import Any

import yaml

from metabolix import __version__
from metabolix.config import DEFAULT_CONFIG, ConfigError, dump_yaml, load_config
from metabolix.pipeline import run_group, run_participant


class Colors:
    """ANSI colors used by the CLI help formatter."""

    HEADER = "\033[95m"
    BLUE = "\033[94m"
    CYAN = "\033[96m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    END = "\033[0m"


class ColoredHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Color section headings and the usage label like the application template."""

    def __init__(self, prog: str, indent_increment: int = 2, max_help_position: int = 36, width: int = 100):
        super().__init__(prog, indent_increment, max_help_position, width)

    def _format_usage(self, usage, actions, groups, prefix):
        if prefix is None:
            prefix = f"{Colors.BOLD}{Colors.GREEN}Usage:{Colors.END} "
        return super()._format_usage(usage, actions, groups, prefix)

    def start_section(self, heading):
        if heading:
            heading = f"{Colors.BOLD}{Colors.CYAN}{heading}{Colors.END}"
        super().start_section(heading)


class ColoredLogFormatter(logging.Formatter):
    """Color runtime log levels on terminals while respecting NO_COLOR."""

    LEVEL_COLORS = {
        logging.DEBUG: Colors.BLUE,
        logging.INFO: Colors.GREEN,
        logging.WARNING: Colors.YELLOW,
        logging.ERROR: Colors.RED,
        logging.CRITICAL: f"{Colors.BOLD}{Colors.RED}",
    }

    def __init__(self, use_color: bool):
        super().__init__("%(asctime)s %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        self.use_color = use_color

    def format(self, record: logging.LogRecord) -> str:
        formatted = super().format(record)
        color = self.LEVEL_COLORS.get(record.levelno)
        if self.use_color and color:
            marker = f" {record.levelname} "
            start = formatted.find(marker)
            if start >= 0:
                end = start + len(marker)
                formatted = f"{formatted[:start]} {color}{record.levelname}{Colors.END} {formatted[end:]}"
        return formatted


def configure_logging(verbose: bool) -> logging.Logger:
    handler = logging.StreamHandler()
    handler.setFormatter(ColoredLogFormatter(use_color=sys.stderr.isatty() and "NO_COLOR" not in os.environ))
    root = logging.getLogger()
    root.handlers.clear()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG if verbose else logging.INFO)
    return logging.getLogger("metabolix")


def build_parser() -> argparse.ArgumentParser:
    description = textwrap.dedent(
        f"""\
        {Colors.BOLD}{Colors.GREEN}METABOLIX v{__version__}
        BIDS-based multi-voxel MRSI application{Colors.END}

        {Colors.BOLD}Description:{Colors.END}
          Process already-converted, image-space MRSI with FSL-MRS. All valid spatial
          voxels are included by default; optional anatomy overlap never restricts fitting.

        {Colors.BOLD}Workflow:{Colors.END}
          1. Discover and validate NIfTI-MRS inputs
          2. Combine tagged receiver coils when needed
          3. Optionally compare residual-water removal branches
          4. Fit all valid spatial voxels and export QC/provenance
          5. Optionally summarize anatomical ROI overlap
        """
    )
    epilog = textwrap.dedent(
        f"""\
        {Colors.BOLD}{Colors.GREEN}EXAMPLES{Colors.END}

        {Colors.BOLD}Dry-run a participant:{Colors.END}
          {Colors.YELLOW}metabolix /data/bids /data/derivatives/metabolix participant \\
            --participant-label 11 --config mrsi.yaml --dry-run{Colors.END}

        {Colors.BOLD}Fit all valid voxels:{Colors.END}
          {Colors.YELLOW}metabolix /data/bids /data/derivatives/metabolix participant \\
            --participant-label 11 --basis /data/basis.BASIS{Colors.END}

        {Colors.BOLD}Compare water-removal branches and report thalamus overlap:{Colors.END}
          {Colors.YELLOW}metabolix /data/bids /data/derivatives/metabolix participant \\
            --water-removal compare --freesurfer-dir /data/derivatives/freesurfer \\
            --roi-labels 10 49 --n-jobs 4{Colors.END}

        {Colors.BOLD}Group aggregation:{Colors.END}
          {Colors.YELLOW}metabolix /data/bids /data/derivatives/metabolix group{Colors.END}

        {Colors.BOLD}{Colors.GREEN}SCIENTIFIC NOTE{Colors.END}
          Successful processing does not establish validated metabolite measurements.
          See README.md and docs/LIMITATIONS.md before interpreting estimates.
        """
    )
    parser = argparse.ArgumentParser(
        prog="metabolix",
        description=description,
        epilog=epilog,
        formatter_class=ColoredHelpFormatter,
        add_help=False,
    )

    required = parser.add_argument_group(f"{Colors.BOLD}Required Arguments{Colors.END}")
    required.add_argument("bids_dir", nargs="?", type=Path, metavar="BIDS_DIR", help="BIDS dataset root (read-only).")
    required.add_argument("output_dir", nargs="?", type=Path, metavar="OUTPUT_DIR", help="Output derivatives directory.")
    required.add_argument("analysis_level", nargs="?", choices=("participant", "group"), metavar="{participant,group}", help="Analysis level.")

    general = parser.add_argument_group(f"{Colors.BOLD}General Options{Colors.END}")
    general.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    general.add_argument("--version", action="version", version=f"metabolix {__version__}")
    general.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging.")
    general.add_argument("-c", "--config", type=Path, metavar="FILE", help="YAML configuration; explicit CLI values override it.")
    general.add_argument("--check-deps", action="store_true", help="Report availability of required FSL-MRS executables.")
    general.add_argument("--export-config", type=Path, metavar="PATH", help="Write the resolved configuration to YAML.")
    general.add_argument("--write-default-config", type=Path, metavar="PATH", help="Write a neutral default YAML configuration.")
    general.add_argument("--write-example-config", type=Path, metavar="PATH", help="Write the illustrative GE PRESS TE144 YAML profile.")

    filters = parser.add_argument_group(
        f"{Colors.BOLD}BIDS Entity Filters{Colors.END}",
        "Select which MRS acquisitions to process.",
    )
    filters.add_argument("-p", "--participant-label", nargs="+", default=None, metavar="LABEL", help="Participant labels, with or without sub- prefix.")
    filters.add_argument("-s", "--session", nargs="+", default=None, metavar="LABEL", help="Session labels.")
    filters.add_argument("--acquisition", nargs="+", default=None, metavar="LABEL", help="Acquisition labels.")
    filters.add_argument("-r", "--run", nargs="+", default=None, metavar="LABEL", help="Run labels.")

    processing = parser.add_argument_group(f"{Colors.BOLD}Processing Options{Colors.END}")
    processing.add_argument("--basis", type=Path, metavar="PATH", help="Explicit compatible FSL-MRS basis.")
    processing.add_argument("--binary-dir", type=Path, metavar="PATH", help="Directory containing FSL-MRS executables.")
    processing.add_argument("--coil-combine", choices=("auto", "on", "off"), help="Coil-combination behavior (default: infer from NIfTI-MRS tags).")
    processing.add_argument("--coil-reference", type=Path, metavar="PATH", help="Optional validated MRSI coil reference.")
    processing.add_argument("--water-removal", choices=("off", "on", "compare"), help="Residual-water handling mode.")
    processing.add_argument("--water-ppm", nargs=2, type=float, metavar=("LOW", "HIGH"), help="Water-removal interval in ppm.")
    processing.add_argument("--ppmlim", nargs=2, type=float, metavar=("LOW", "HIGH"), help="Fitting interval in ppm.")
    processing.add_argument("--ignore", nargs="+", metavar="COMPONENT", help="Basis components to exclude from fitting.")
    processing.add_argument("--combine", nargs=2, action="append", metavar=("A", "B"), help="Combine component estimates for reporting; repeatable.")
    processing.add_argument("--internal-ref", nargs="+", metavar="COMPONENT", help="Components used as the internal reference.")
    processing.add_argument("--baseline", metavar="MODEL", help="FSL-MRS baseline model.")
    processing.add_argument("--algorithm", choices=("Newton", "MH"), help="FSL-MRS fitting algorithm.")
    processing.add_argument("--fit-mask", type=Path, metavar="PATH", help="Optional spatial fit mask matching the MRS grid.")
    processing.add_argument("--n-jobs", type=int, metavar="N", help="Bound FSL-MRS local parallel workers.")

    anatomy = parser.add_argument_group(f"{Colors.BOLD}Anatomical ROI Options{Colors.END}")
    anatomy.add_argument("--freesurfer-dir", type=Path, metavar="PATH", help="FreeSurfer derivatives root containing aseg.mgz.")
    anatomy.add_argument("--roi-map", type=Path, metavar="PATH", help="Integer label map or binary ROI mask.")
    anatomy.add_argument("--roi-labels", nargs="+", type=int, metavar="INTEGER", help="Labels to summarize (10 49 names left/right thalamus).")
    anatomy.add_argument("--roi-step-mm", type=float, metavar="MM", help="Maximum ROI sampling spacing (<= 0.5 mm).")
    anatomy.add_argument("--minimum-overlap-percent", type=float, metavar="PCT", help="Reporting-only overlap threshold.")

    reporting = parser.add_argument_group(f"{Colors.BOLD}Reporting and Execution Options{Colors.END}")
    reporting.add_argument("--report-metabolites", nargs="+", metavar="NAME", help="Metabolites to emphasize in report summaries; does not prune fit.")
    reporting.add_argument("--report", action=argparse.BooleanOptionalAction, default=None, help="Enable/disable the FSL-MRS report.")
    reporting.add_argument("--resume", action=argparse.BooleanOptionalAction, default=None, help="Resume only when run provenance matches.")
    reporting.add_argument("--overwrite", action=argparse.BooleanOptionalAction, default=None, help="Explicitly replace an existing output run.")
    reporting.add_argument("--dry-run", action="store_true", help="Validate inputs and print planned commands without running FSL-MRS.")
    return parser


def _overrides(args: argparse.Namespace) -> dict[str, Any]:
    sections: dict[str, dict[str, Any]] = {}

    def put(section: str, key: str, value: Any, transform=lambda item: item) -> None:
        if value is not None:
            sections.setdefault(section, {})[key] = transform(value)

    for arg, key in (("participant_label", "participant_label"), ("session", "session"), ("acquisition", "acquisition"), ("run", "run")):
        put("selection", key, getattr(args, arg))
    put("basis", "path", args.basis, lambda value: str(value.resolve()))
    put("basis", "binary_dir", args.binary_dir, lambda value: str(value.resolve()))
    put("processing", "water_removal", args.water_removal)
    put("processing", "coil_reference", args.coil_reference, lambda value: str(value.resolve()))
    put("processing", "water_ppm", args.water_ppm)
    put("processing", "coil_combine", args.coil_combine)
    put("fit", "ppmlim", args.ppmlim)
    put("fit", "ignore", args.ignore)
    put("fit", "combine", args.combine)
    put("fit", "internal_reference", args.internal_ref)
    put("fit", "baseline", args.baseline)
    put("fit", "algorithm", args.algorithm)
    put("fit", "mask", args.fit_mask, lambda value: str(value.resolve()))
    put("roi", "freesurfer_dir", args.freesurfer_dir, lambda value: str(value.resolve()))
    put("roi", "map", args.roi_map, lambda value: str(value.resolve()))
    if args.roi_labels is not None:
        if args.roi_labels == [10, 49]:
            label_map = {"left_thalamus": 10, "right_thalamus": 49}
        else:
            label_map = {f"label_{label}": label for label in args.roi_labels}
        put("roi", "labels", label_map)
    put("roi", "sampling_step_mm", args.roi_step_mm)
    put("roi", "minimum_overlap_percent", args.minimum_overlap_percent)
    put("execution", "n_jobs", args.n_jobs)
    put("execution", "resume", args.resume)
    put("execution", "overwrite", args.overwrite)
    put("report", "enabled", args.report)
    put("report", "metabolites", args.report_metabolites)
    return sections


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logger = configure_logging(args.verbose)
    try:
        if args.write_default_config:
            dump_yaml(DEFAULT_CONFIG, args.write_default_config)
            logger.info("Wrote neutral config to %s", args.write_default_config)
            if not any((args.write_example_config, args.check_deps, args.bids_dir)):
                return 0
        if args.write_example_config:
            from metabolix.examples import example_config

            dump_yaml(example_config(), args.write_example_config)
            logger.info("Wrote GE PRESS TE144 example config to %s", args.write_example_config)
            if not any((args.check_deps, args.bids_dir)):
                return 0
        config = load_config(args.config, _overrides(args))
        if args.export_config:
            dump_yaml(config, args.export_config)
        if args.check_deps:
            from metabolix.execution import dependency_report

            report = dependency_report(config)
            for name, location in report.items():
                print(f"{name}: {location or 'NOT FOUND'}")
            return 0 if all(report.values()) else 2
        if not (args.bids_dir and args.output_dir and args.analysis_level):
            parser.error("BIDS_DIR OUTPUT_DIR and analysis_level are required unless generating a config.")
        args.bids_dir = args.bids_dir.resolve()
        args.output_dir = args.output_dir.resolve()
        if not args.bids_dir.is_dir():
            raise ConfigError(f"BIDS input directory does not exist: {args.bids_dir}")
        if args.analysis_level == "participant":
            return run_participant(args.bids_dir, args.output_dir, config, args.dry_run, logger)
        return run_group(args.output_dir, config, args.dry_run, logger)
    except (ConfigError, OSError, RuntimeError, ValueError) as exc:
        logger.error("%s", exc)
        return 2


if __name__ == "__main__":
    sys.exit(main())