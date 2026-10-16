# Browser aligner builds

The Pages build compiles MAFFT `tbfast` and `dvtditr` (Biowasm package
7.520) and MUSCLE (Biowasm package 5.1.0) with Emscripten 2.0.25,
`ALLOW_MEMORY_GROWTH=1` and `MAXIMUM_MEMORY=4294967296`. This changes the
wasm32 heap ceiling from the default 2 GiB to 4 GiB per module. It does
not allocate 4 GiB immediately, change the alignment strategy, or promise
that a browser can supply that much memory. The page, Pyodide, files and
other modules also consume memory. A workload can still exceed the limit.

`scripts/build_browser_aligners.py` pins the Biowasm repository and both
source revisions, applies its existing tool patches and adds the larger
heap ceiling. The generated JavaScript, WebAssembly, licenses, patches and
build manifest are published under `site/aligners/`. They are build
artifacts rather than committed binaries. Aioli loads these same-origin
programs via `urlPrefix`; it does not silently fall back to CDN aligner
binaries with the old ceiling. Aioli itself and its small base filesystem
module continue to load from Biowasm.

The package version is retained for provenance. In particular, the MUSCLE
package's own banner can report a different program version. Its RAM banner
is a fixed value in the upstream Biowasm compatibility patch; it is not
the configured heap maximum or the host computer's available RAM.

## Build and verify

Run from the repository root with Docker available:

```bash
docker run --rm -e HOME=/tmp \
  -v "$PWD:/src" -w /src emscripten/emsdk:2.0.25 \
  python3 scripts/build_browser_aligners.py
node site/aligner-build-test.mjs
node site/build-engine.mjs
node site/smoke-test.mjs
```

Alternatively activate Emscripten 2.0.25 and run the Python build script
directly. Use a new `--work-dir` for each rebuild; the script refuses to
overwrite an existing source checkout. Engine bundling also requires the
usual Python/Biopython dependencies. Serve `site/` over HTTP for previews.

CI verifies all three wasm modules' declared maxima and SHA256 checksums,
then uses Chromium to allocate and access memory beyond 2 GiB in a matching
probe build. It also runs two real alignments with each aligner through
the production worker, checking taxon count, aligned lengths and ungapped
sequence preservation. Playwright is only a test dependency. The probe is
for build validation and is not automatically run in users' browser sessions.

MAFFT first-stage and refinement failures stop the workflow immediately.
A failed `tbfast` does not trigger `dvtditr`, and an existing partial `pre`
file cannot make a reported refinement failure appear successful. Memory
messages name the failing aligner, preserve its diagnostic and report an
explicit byte limit when supplied by the runtime. Only MUSCLE failures
suggest trying MAFFT or explain the MUSCLE RAM banner. Aligner selection
never changes automatically.

This is a 4 GiB wasm32 build, not a Memory64 port. It does not guarantee that
the previously failing five-sequence, approximately 10,146-residue dataset
will fit: that dataset needs a separate run to establish its peak demand.

Custom preview servers must serve the `.js` and `.wasm` files with appropriate MIME types and `Access-Control-Allow-Origin: *`, because Aioli can load its nested worker through a data URL. The browser test's local server sets these headers. Tests can optionally mirror the three unchanged CDN runtime assets through `PHYLOPREP_TEST_VENDOR_DIR` when the test browser cannot use a host network proxy; production aligner code is not rewritten.
