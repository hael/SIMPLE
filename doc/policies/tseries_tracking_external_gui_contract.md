# SINGLE–Nanotation tracking interface

Scope: the full-field display images and project metadata written by
`single_exec prg=tseries_prep4tracking`, their annotation in the external
[Nanotation GUI](https://github.com/cyrilreboul/nanotation), and the mapping of
returned annotations to original frames.

**Current status:** image preparation is available in SINGLE. The historical
annotation-ingestion application, `tseries_extractor`, remains retired. The
outgoing handoff and the return-file format below do not imply that a complete
GUI-to-SINGLE extraction workflow is currently executable.

The GUI behavior described here was checked against Nanotation revision
[`baf6092f8df96ca221e77200765683766684ce24`](https://github.com/cyrilreboul/nanotation/commit/baf6092f8df96ca221e77200765683766684ce24).
Recheck the upstream exporter and coordinate conventions when changing either
side of this interface.

## Workflow and ownership

```text
Original frames + SIMPLE project
    -> tseries_prep4tracking
    -> preparation project copy + filtered/*.mrc
    -> Nanotation opens the filtered/ directory
    -> annotation text files with image paths and selected coordinates
    -> SINGLE frame mapping and trajectory extraction (historical consumer)
```

SINGLE owns the original frame identities, project rows, display-image generation
and mapping between display groups and original frames. Nanotation owns the
interactive annotations and their text export. It is maintained as a separate
repository and application.

Nanotation opens MRC files from a directory. It does **not** read the `.simple`
project or its `track_fname`, `fromt` and `tot` fields. The exported image filename
is the join key back to SINGLE's project metadata. These fields are therefore an
external-workflow contract even when no current SIMPLE routine reads them.
Removing or renaming them requires reviewing the Nanotation return handoff;
absence of an in-repository reader is not evidence that the metadata is disposable.

## Preparing and opening the display images

```sh
single_exec prg=tseries_prep4tracking projfile=/path/to/project.simple nthr=8
```

The command creates its own execution directory and project copy. Keep that copy:
it contains the mapping for the generated display images. The original project
and original frame files are unchanged. The frames are read from the project's
`mic.frame` records, in project-row order.

Inside the preparation directory, SINGLE writes:

| Output | Purpose |
| --- | --- |
| `filtered/frame_<index>.mrc` | Processed, potentially downscaled full-field images for annotation |
| The copied `.simple` project | `mic` records linking those images to original frame ranges |
| `global_avgimg.mrc` | Native-grid average of every original frame, for diagnostics |
| `outliers.txt` | One-based native-grid hot/dead-pixel coordinates, for diagnostics |

Open the `filtered/` directory in Nanotation:

```sh
nanotation /path/to/preparation_run/filtered
```

Use the GUI's own [README and installation instructions](https://github.com/cyrilreboul/nanotation#readme).
The GUI has its own Python environment and dependencies. Open only the intended
display-image directory; its loaded-image ordinal is separate from SINGLE's
original frame index.

## Project metadata schema

All fields below belong to the **`mic` segment**, on the row for the first original
frame of a display group. They are not particle-orientation slots. For a project
with `N` frames, group-start rows are `s = 1, 11, 21, ... <= N`.

| Field | Type | Writer | Meaning |
| --- | --- | --- | --- |
| `track_fname` | String: absolute canonical MRC pathname | `tseries_prep4tracking` | Identity of the generated display image. The returned annotation's `filename` must match this value |
| `fromt` | Integer, one-based | `tseries_prep4tracking` | First original frame associated with this group: `s` |
| `tot` | Integer, one-based and inclusive | `tseries_prep4tracking` | Last original frame associated with this group: `min(s + 9, N)` |

The grouping range and final display-summation window differ. Each image sums
up to five **denoised frame estimates** beginning at `s`, within its batch;
its metadata covers up to ten original frames. For example, `frame_001.mrc`
sums denoised estimates 1–5 and has `fromt=1`, `tot=10`. A final group beginning
at 101 in a 102-frame project sums estimates 101–102 and has `fromt=101`,
`tot=102`. Each estimate can borrow information from elsewhere in its FLEX
batch, subject to the temporal/parity controls; the summation window does not
identify all contributing raw pixels.

Images are prepared in batches of at most 100, with background subtraction, FLEX
kernel denoising, short-frame summation, nonlocal-means filtering, bandpass and
normalization. Their pixels are display products; the historical extraction
reader used native-resolution `mic.intg` images when present and otherwise
`mic.frame`. The GUI's selected positions must be interpreted on the display grid.

## Returned annotations and coordinate conventions

Use Nanotation's **Export track** action to write the text interchange (default
name `annotations.txt`). Saving a GUI session writes a separate JSON document;
it is not the SINGLE coordinate interchange.

The text exporter writes one space-separated `key=value` record per loaded
display image between the first and last annotation checkpoints. Positions can
be interpolated or smoothed, including at checkpoint images; the rows need not
be raw clicks. The checked GUI revision defaults to smoothness `0.15`, uses
linear interpolation for two distinct checkpoints and a spline when at least
three are available. An export represents one connected track and has no
particle-ID column. The record includes the resolved image pathname, position,
loaded-image ordinal and displayed image dimensions. For example:

```text
filename=/data/run/1_tseries_prep4tracking/filtered/frame_001.mrc x=250.5 y=301.2 index=1 xdim=1024 ydim=768
```

| Field | Meaning |
| --- | --- |
| `filename` | Absolute resolved pathname of the annotated display image; join to `mic.track_fname` |
| `x`, `y` | Zero-based selected centre coordinates on the displayed image, with a bottom-left origin |
| `index` | One-based ordinal in Nanotation's loaded sequence; not the original frame number |
| `xdim`, `ydim` | Displayed image width and height in pixels; not the native frame dimensions |

For files named `frame_001.mrc`, `frame_011.mrc`, `frame_021.mrc`, the GUI ordinals
are 1, 2, 3 while the corresponding original frame starts are 1, 11, 21. The GUI
exports no `fromt`, `tot` or native pixel-size metadata;
the prepared project supplies those associations.

Retain the exported path and dimensions. Do not substitute `index` for `fromt`
or use the display coordinates directly as native extraction-box corners.
For display dimensions `(Dx,Dy)` and native dimensions `(Nx,Ny)`, the historical
SINGLE reader used:

```text
xscale = Dx / Nx                     yscale = Dy / Ny
native centre x = exported x/xscale  native centre y = exported y/yscale
native box corner = native centre - box_extract/2
```

It rounded the resulting native box corners to integer pixels and interpolated
between successive selected positions over the project frame ranges. The final
selection's position was held across its remaining range. No additional Y-axis
flip was applied by that reader; the export already used the bottom-left
coordinate convention. Nanotation's viewer orientation makes Y increase upward,
but the exporter applies no numerical Y flip or half-pixel offset and writes
image-grid pixel values rather than physical Angstrom coordinates. Coordinates
are not clamped to the image bounds by the checked exporter; a return consumer
must validate extraction bounds.

## Path identity and file lifetime

Keep the preparation project copy, display images, original frames and annotation
files together as the associated dataset. The historical reader compared returned
`filename` and project `track_fname` **as exact strings**, in increasing frame order.
SINGLE canonicalizes existing display-image paths, and Nanotation resolves them
when exporting, so an ordinary symlink alias to the same file resolves to the
same identity. Moving or renaming the files, or retargeting an alias to another
location, can change that identity even when the pixels are identical.

There is no automatic remapping of these fields by this application. A relocation
or future importer must update or validate both sides explicitly. The historical
SINGLE text parser separates fields on spaces and has no quoted-path support;
use paths without whitespace for that interchange format.

## SINGLE re-entry status

The historical return command was:

```text
single_exec prg=tseries_extractor projfile=<prepared-project-copy> infile=<export.txt> box_extract=<native-pixel-box>
```

This is a description of the historical interface, **not a command available in
the current executable**. Its implementation at `553ba21c4` matched `filename` to
`track_fname`, mapped through `fromt`/`tot`, rescaled the positions and extracted
from native-resolution `mic.intg` images when present, falling back to `mic.frame`.
It was removed in `703a7e825` and has not been restored with the image-preparation
application.

`tseries_track_particles` is the automated tracker taking an EMAN box file;
`trajectory_import` imports an already extracted image stack. Neither is a
replacement reader for Nanotation's exported coordinate records. A restored or
replacement return consumer must explicitly implement and validate the mapping
above before the round trip is described as supported.

## Compatibility and verification

Preserve key names, absolute path identity, one-based inclusive original-frame
ranges, displayed-grid coordinate meaning and the distinction between GUI
ordinals and original frame indices. Review changes to these conventions with
the Nanotation maintainer. Do not change the interface merely because a field
has no active reader within this repository.

SINGLE's native preparation tests validate project persistence, group ranges,
display-image headers, proportional rectangular scaling and source preservation.
They do not launch Nanotation or validate a currently available return consumer.
An end-to-end compatibility check needs an export from the identified GUI
revision and a working SINGLE annotation reader.

## Implementation references

- [Nanotation repository](https://github.com/cyrilreboul/nanotation) and its
  [README](https://github.com/cyrilreboul/nanotation#readme).
- Checked GUI sources: [launch and coordinates](https://github.com/cyrilreboul/nanotation/blob/baf6092f8df96ca221e77200765683766684ce24/README.md#L33-L67),
  [folder reader](https://github.com/cyrilreboul/nanotation/blob/baf6092f8df96ca221e77200765683766684ce24/src/nanotation/mrc_io.py#L65-L100),
  [text exporter](https://github.com/cyrilreboul/nanotation/blob/baf6092f8df96ca221e77200765683766684ce24/src/nanotation/annotations.py#L130-L187),
  [viewer orientation](https://github.com/cyrilreboul/nanotation/blob/baf6092f8df96ca221e77200765683766684ce24/src/nanotation/napari_ui.py#L7-L23)
  and [export action](https://github.com/cyrilreboul/nanotation/blob/baf6092f8df96ca221e77200765683766684ce24/src/nanotation/widgets.py#L637-L666).
- [SINGLE producer](../../src/main/strategies/parallelization/simple_tseries_prep_strategy.f90),
  [CLI metadata](../../src/main/ui/single/single_ui_tseries.f90) and
  [preparation tests](../../src/main/commanders/single/simple_tseries_prep4tracking_tester.f90).
