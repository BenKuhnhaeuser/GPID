# Conda recipe for GPID

This recipe builds GPID as a script-based conda package. The package installs
the user-facing `gpid` command into `${PREFIX}/bin` and stores workflow scripts
and templates under `${PREFIX}/share/gpid`.

Build locally with:

```bash
conda build conda-recipe -c conda-forge -c bioconda
```

The recipe currently uses `source: path: ..` for local builds from this working
tree. For a public feedstock, create a versioned source archive, update
`source` to the archive URL, and add the SHA-256 checksum.
