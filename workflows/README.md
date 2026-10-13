# Workflow definitions

Snakemake/Nextflow entrypoints that chain the stages, for running the pipeline
non-interactively on a cluster.

Each stage is independently runnable via the `phylodonor` CLI, which is the supported
interface. Workflows here are a convenience layer over it, not a parallel implementation —
they shell out to the same CLI so there is one code path.
