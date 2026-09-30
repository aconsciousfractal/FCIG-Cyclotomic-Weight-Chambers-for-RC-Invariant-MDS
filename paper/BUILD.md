# Building the article

From the repository root:

```bash
python -X utf8 -B scripts/build_paper.py
```

The build requires Tectonic 0.16.9. Set `TECTONIC` to an explicit executable
if needed. The script compiles twice in clean temporary directories, requires
byte-identical output inside that toolchain, rejects missing citations,
undefined references, and overfull boxes, then writes the title-named PDF and
`results/paper_build.json`.

For a non-writing, byte-exact comparison with the committed PDF and its full
build receipt:

```bash
python -X utf8 -B scripts/build_paper.py --verify-rebuild
```

This route requires the exact `Tectonic 0.16.9` version string, rebuilds twice,
compares the rebuilt bytes with the committed PDF, and recomputes every receipt
field, including the PDF SHA-256, source-tree SHA-256, theorem labels, page
count, byte count, toolchain, and source-date epoch.

Versions 3.0.0, 3.1.0 and 3.2.0 use the available Tectonic 0.16.9 toolchain. This explicit toolchain change supersedes the v2 PDF build environment; byte identity is required between clean v3 builds, not between versions. The source is a multi-file LaTeX project.
