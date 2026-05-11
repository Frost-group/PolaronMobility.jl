# Documentation

The documentation is built from the Markdown files in `docs/src` plus exported
Julia docstrings through Documenter.jl.

## Local build

```bash
julia --project=docs docs/make.jl
```

For a local live preview of the generated static site, serve `docs/build` with
a simple static-file server after the Documenter build completes.

The source pages of interest for the lattice theory are:

- `docs/src/scientific_discussion.md`
- `docs/src/lattice_transport.md`
- `docs/src/examples.md`
- `docs/src/functions.md`

The active lattice transport story in those pages is the CTMC first-return
current-blip / exact-sideband formulation for Holstein, Peierls, and
Holstein-Peierls models. The repository README summarizes the same public
theory and API at a higher level.
