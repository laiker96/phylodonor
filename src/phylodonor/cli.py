"""Command-line entrypoint.

    phylodonor <stage> <step> --config configs/my_target.yaml

Every stage takes the same ``--config``, validates it before doing any work,
and writes under ``results/<project.name>/<stage>/``. Stages are deliberately
separate commands rather than one monolithic run: the expensive steps
(alignment, tree inference, simulation) are ones you re-run selectively.

Use ``--dry-run`` to validate a config and print the resolved inputs and
output paths without executing anything. That is the fastest way to find out
that a database path is wrong.
"""

from __future__ import annotations

import argparse
import sys

from phylodonor import __version__, config as config_mod

STAGES: dict[str, list[str]] = {
    "homologs": ["search", "annotate", "classify"],
    "phylogeny": ["align", "qc", "tree"],
    "thermal": ["annotate", "clades"],
    "donors": ["map", "score", "rank"],
    "holo": ["build", "minimize", "qc"],
    "grafting": ["enumerate", "score", "rank"],
    "dynamics": ["prepare", "run", "analyze"],
    "scoring": ["conservation", "solubility", "aggregation"],
    "viz": ["tree", "logos", "summary"],
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="phylodonor",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"phylodonor {__version__}")
    subparsers = parser.add_subparsers(dest="stage", metavar="stage", required=True)

    for stage, steps in STAGES.items():
        sp = subparsers.add_parser(stage, help=f"{stage} stage ({', '.join(steps)})")
        sp.add_argument("step", choices=steps, help="step to run")
        sp.add_argument("-c", "--config", required=True, help="target YAML config")
        sp.add_argument("--dry-run", action="store_true",
                        help="validate config and report paths without running")
        sp.add_argument("-v", "--verbose", action="store_true")

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    try:
        cfg = config_mod.load(args.config)
        cfg.check_stage(args.stage)
    except config_mod.ConfigError as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return 2

    out = cfg.results_dir(args.stage)
    if args.dry_run:
        print(f"config   {cfg.path}")
        print(f"project  {cfg.name}")
        print(f"stage    {args.stage}:{args.step}")
        print(f"outputs  {out}")
        return 0

    # Stage dispatch is wired up as each stage is ported. Until then, fail
    # loudly rather than appearing to succeed.
    print(f"not yet implemented: {args.stage}:{args.step}", file=sys.stderr)
    print("See docs/roadmap.md for the migration order.", file=sys.stderr)
    return 3


if __name__ == "__main__":
    sys.exit(main())
