from pathlib import Path

from django.template.loader import render_to_string
from django.test import SimpleTestCase


class BatchDetailTemplateTests(SimpleTestCase):
    @staticmethod
    def _read_template(relative_path):
        base_dir = Path(__file__).resolve().parents[1]
        return (base_dir / "templates" / relative_path).read_text(encoding="utf-8")

    @staticmethod
    def _read_static(relative_path):
        base_dir = Path(__file__).resolve().parents[1]
        return (base_dir / "static" / relative_path).read_text(encoding="utf-8")

    @staticmethod
    def _preview(
        path,
        name,
        kind,
        crop,
        hidden_by_default=False,
        visibility_group=None,
        width=None,
        height=None,
    ):
        label = kind.replace("_", " ")
        preview = {
            "path": path,
            "name": name,
            "kind": kind,
            "crop": crop,
            "title": f"{name} — {label}",
            "alt": f"{name} {label}",
        }
        if hidden_by_default:
            preview["hidden_by_default"] = True
        if visibility_group is not None:
            preview["visibility_group"] = visibility_group
        if width is not None and height is not None:
            preview["width"] = width
            preview["height"] = height
        return preview

    def test_batch_detail_template_uses_mobile_first_layout(self):
        batch_view = self._read_template("nice_classic/batchview.html")

        self.assertIn("flex flex-wrap sm:flex-nowrap", batch_view)
        self.assertIn('aria-label="back to workspace"', batch_view)
        self.assertIn('aria-label="stop batch job"', batch_view)
        self.assertIn('aria-label="delete batch job"', batch_view)
        self.assertIn("size-11", batch_view)
        self.assertIn("grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3", batch_view)
        self.assertIn('data-collapse-key="artifacts" class="lg:col-span-full', batch_view)
        self.assertIn("flex overflow-x-auto rounded-t-lg", batch_view)
        self.assertIn("grid grid-cols-2 sm:grid-cols-3 gap-2", batch_view)
        self.assertIn('id="batch_artifact_images" data-output-grid', batch_view)
        self.assertIn("h-[220px] sm:h-[260px]", batch_view)
        self.assertIn("minmax(0,1fr)", batch_view)
        self.assertIn("touchmove", batch_view)
        self.assertIn("capture: true", batch_view)

    def test_batch_detail_reuses_stream_pills_and_composite_crop_pattern(self):
        batch_view = self._read_template("nice_classic/batchview.html")
        preview = self._read_template("nice_classic/includes/_batch_artifact_preview.html")

        self.assertIn("nice_stream/includes/_pill.html", batch_view)
        self.assertIn("nice_classic/includes/_batch_artifact_preview.html", batch_view)
        self.assertIn("absolute inset-0 w-full h-full object-cover", preview)
        self.assertIn("object-[100%]", preview)
        self.assertIn('style="object-position:0%"', preview)

    def test_batch_detail_reuses_stream_class_selector_component(self):
        batch_view = self._read_template("nice_classic/batchview.html")
        selector = self._read_template("nice_stream/includes/_class_selector.html")
        display_levels = self._read_template(
            "nice_classic/includes/_batch_display_levels.html"
        )
        selector_script = self._read_static("nice_lite/class_selector.js")
        slider_script = self._read_static("nice_lite/pick_micrograph_slider.js")
        output_grid_script = self._read_static("nice_lite/output_grid.js")
        rendered = render_to_string(
            "nice_classic/batchview.html",
            {
                "jobid": 7,
                "status": "finished",
                "result_project": "/project/workspace.simple",
                "project_summary_available": False,
                "class_selector_available": True,
                "class_selector_requested": True,
                "batch_class_selector": {
                    "classes": ({
                        "class_id": 1,
                        "stack_index": 1,
                        "population": 25,
                        "resolution": 4.5,
                    },),
                    "class_count": 1,
                    "stack_name": "classes.mrcs",
                    "width": 128,
                    "height": 128,
                    "initial_selected_class_ids": (1,),
                    "browser_data": {
                        "classes": [{"class_id": 1, "stack_index": 1}],
                        "initial_selected_class_ids": [1],
                        "selection_storage_key": "selection-key",
                    },
                },
                "class_selector_replaces_artifact_previews": True,
                "arguments": [],
                "logs": [],
                "artifact_counts": [],
                "artifact_images": [{
                    "name": "redundant-cavgs.jpg",
                    "previews": [self._preview(
                        "/project/redundant-cavgs.jpg",
                        "redundant-cavgs.jpg",
                        "image",
                        "full",
                    )],
                }],
                "auto_refresh": False,
            },
        )

        self.assertNotIn("data-class-selector-action", batch_view)
        self.assertNotIn("select 2D classes", rendered)
        self.assertNotIn("2D class selector", rendered)
        self.assertNotIn(
            "Review and update the class selection directly from this output.",
            rendered,
        )
        self.assertIn("nice_stream/includes/_class_selector.html", batch_view)
        self.assertLess(
            rendered.index('id="batch_artifacts_panel"'),
            rendered.index('id="batch_class_selector"'),
        )
        self.assertNotIn("redundant-cavgs.jpg", rendered)
        self.assertIn("data-class-selector-grid", selector)
        self.assertIn("data-output-grid", selector)
        self.assertIn('data-output-grid-size-property="--class-tile-size"', selector)
        self.assertIn("data-output-grid-tile", selector)
        self.assertNotIn('id="batch_class_tile_size"', selector)
        self.assertIn('id="batch_output_tile_size" type="range"', rendered)
        self.assertIn('data-output-grid-target-id="batch_artifacts_panel"', rendered)
        self.assertIn("data-class-selector-select-all", selector)
        self.assertIn(">invert selection</button>", rendered)
        self.assertIn(">reset selection</button>", rendered)
        self.assertIn("batch_class_deselection_export", selector)
        self.assertIn('<option value="resolution" selected>resolution</option>', selector)
        self.assertIn("nice_classic/includes/_batch_display_levels.html", selector)
        self.assertIn("data-class-selector-toolbar", selector)
        self.assertIn("data-class-selector-control-row", selector)
        self.assertIn("data-class-selector-sort-field", selector)
        self.assertIn("data-class-selector-bulk-actions", selector)
        self.assertIn("embedded=True", selector)
        self.assertIn(
            "[data-class-selector-sort-field] {\n"
            "                display: flex;\n"
            "                align-items: center;\n"
            "                gap: 0.5rem;",
            batch_view,
        )
        self.assertLess(
            selector.index("data-class-selector-toolbar"),
            selector.index("data-class-selector-sort-field"),
        )
        self.assertLess(
            selector.index("data-class-selector-sort-field"),
            selector.index("nice_classic/includes/_batch_display_levels.html"),
        )
        self.assertLess(
            rendered.index("<label data-class-selector-sort-field"),
            rendered.index(
                "<div data-pick-display-controls data-display-levels "
                "data-display-levels-embedded"
            ),
        )
        self.assertLess(
            rendered.index(
                "<div data-pick-display-controls data-display-levels "
                "data-display-levels-embedded"
            ),
            rendered.index("<div data-class-selector-bulk-actions"),
        )
        self.assertIn('id="batch_class_display_min" type="range"', rendered)
        self.assertIn('id="batch_class_display_max" type="range"', rendered)
        self.assertIn('id="batch_class_display_brightness" type="range"', rendered)
        self.assertIn('id="batch_class_display_contrast" type="range"', rendered)
        self.assertIn('data-pick-display-target-id="batch_class_selector_grid"', rendered)
        self.assertIn('data-pick-display-image-selector="[data-class-selector-tile] img"', rendered)
        self.assertIn('id="batch_pick_overlay_toggle"', rendered)
        self.assertIn('data-pick-overlay-mode="circle"', rendered)
        self.assertIn('data-pick-overlay-modes="circle,boxes"', rendered)
        self.assertIn('data-pick-overlay-name="class-average"', rendered)
        self.assertIn('aria-controls="batch_class_selector_grid"', rendered)
        self.assertIn('aria-label="circle overlay; show box"', rendered)
        self.assertIn('data-pick-overlay-toggle-label>circle', rendered)
        self.assertIn('id="batch_pick_overlay_size" type="range"', rendered)
        self.assertLess(
            rendered.index('id="batch_pick_overlay_size"'),
            rendered.index('id="batch_pick_overlay_toggle"'),
        )
        class_size_input = rendered.split(
            'id="batch_pick_overlay_size" type="range"', 1,
        )[1].split(">", 1)[0]
        self.assertNotIn("disabled", class_size_input)
        self.assertIn('id="batch_pick_overlay_color_control"', rendered)
        self.assertIn('data-pick-overlay-color="#ff3b30"', rendered)
        self.assertIn("data-class-average-overlay-preview", rendered)
        self.assertIn("data-class-average-overlay", rendered)
        self.assertIn('data-source-width="128"', rendered)
        self.assertIn('data-source-height="128"', rendered)
        self.assertIn("drawClassAverageOverlay", slider_script)
        self.assertIn("initializeClassAverageOverlay", slider_script)
        self.assertIn("scaledOverlaySize / 2", slider_script)
        self.assertIn("classAverageResizeObserver", slider_script)
        self.assertIn('mode === "boxes" ? "box" : mode', slider_script)
        self.assertEqual(display_levels.count("accent-color:var(--color-streamaccent)"), 4)
        self.assertEqual(display_levels.count("data-display-level-field"), 4)
        self.assertEqual(display_levels.count("data-display-level-caption"), 4)
        self.assertEqual(display_levels.count("data-display-level-name"), 4)
        self.assertIn("data-display-levels", display_levels)
        self.assertIn("data-display-levels-embedded", rendered)
        self.assertIn("data-display-level-header", display_levels)
        self.assertNotIn("data-display-level-title", display_levels)
        self.assertNotIn("data-display-level-help", display_levels)
        self.assertNotIn(">display levels<", rendered)
        self.assertNotIn(
            "brightness/contrast and min/max are synchronized views of the "
            "displayed 8-bit level window",
            rendered,
        )
        self.assertIn("@media (min-width: 64rem)", batch_view)
        self.assertIn(
            "grid-template-columns: minmax(0, 1fr) max-content",
            batch_view,
        )
        self.assertIn(
            "grid-template-columns: 5rem minmax(0, 1fr) 3rem",
            batch_view,
        )
        self.assertIn(
            '[data-display-level-field] > input[type="range"] {\n'
            "                grid-column: 2;\n"
            "                grid-row: 1;",
            batch_view,
        )
        self.assertIn(
            "[data-display-level-caption] {\n"
            "                display: contents;",
            batch_view,
        )
        self.assertIn(
            "[data-display-level-name] {\n"
            "                grid-column: 1;\n"
            "                grid-row: 1;\n"
            "                justify-self: end;",
            batch_view,
        )
        self.assertIn(
            "[data-display-level-field] [data-pick-display-output] {\n"
            "                grid-column: 3;\n"
            "                grid-row: 1;\n"
            "                text-align: left;",
            batch_view,
        )
        self.assertIn("text-align: right", batch_view)
        self.assertIn(
            "grid-template-columns: max-content minmax(0, 1fr) max-content",
            batch_view,
        )
        self.assertIn("data-display-level-grid", display_levels)
        self.assertIn(
            "grid-template-columns: repeat(2, minmax(0, 1fr))",
            batch_view,
        )
        self.assertIn("@media (min-width: 72rem)", batch_view)
        self.assertIn(
            "grid-template-columns: repeat(4, minmax(0, 1fr))",
            batch_view,
        )
        self.assertLess(
            batch_view.index("@media (min-width: 64rem)"),
            batch_view.index("@media (min-width: 72rem)"),
        )
        self.assertLess(
            display_levels.index("_display_brightness"),
            display_levels.index("_display_contrast"),
        )
        self.assertLess(
            display_levels.index("_display_contrast"),
            display_levels.index("_display_min"),
        )
        self.assertLess(
            display_levels.index("_display_min"),
            display_levels.index("_display_max"),
        )
        self.assertIn("128 × 128 px", rendered)
        self.assertIn("applyVisualRange", selector_script)
        self.assertIn('let sortKey = sortControl.value || "resolution"', selector_script)
        self.assertIn("localStorage", selector_script)
        self.assertNotIn("data-class-selector-brightness", selector_script)
        self.assertIn("pickDisplayImageSelector", slider_script)
        self.assertIn("target.querySelectorAll(imageSelector)", slider_script)
        self.assertNotIn("wheelAdjustedTileSize", selector_script)
        self.assertIn("wheelAdjustedTileSize", output_grid_script)
        self.assertIn("event.ctrlKey", output_grid_script)
        self.assertIn("event.metaKey", output_grid_script)
        self.assertIn("event.preventDefault()", output_grid_script)
        self.assertIn("{ passive: false }", output_grid_script)
        self.assertIn('target.querySelectorAll("[data-output-grid]")', output_grid_script)
        self.assertIn("output_grid.js?v=2", rendered)
        self.assertIn("pick_micrograph_slider.js?v=14", rendered)
        self.assertIn("class_selector.js?v=4", rendered)

    def test_batch_detail_ctf_artifacts_add_source_micrographs_automatically(self):
        batch_view = self._read_template("nice_classic/batchview.html")
        diagnostic_path = "/project/workspace/2_ctf_estimate/movie_ctf_estimate_diag.jpg"
        micrograph_path = "/project/workspace/1_motion_correct/movie_thumb.jpg"
        context = {
            "jobid": 7,
            "artifact_counts": [{"extension": "JPG", "count": 1}],
            "artifact_images": [{
                "name": "movie_ctf_estimate_diag.jpg",
                "path": diagnostic_path,
                "previews": [
                    self._preview(
                        diagnostic_path,
                        "movie_ctf_estimate_diag.jpg",
                        "image",
                        "full",
                        width=800,
                        height=600,
                    ),
                    self._preview(
                        micrograph_path,
                        "movie_thumb.jpg",
                        "micrograph",
                        "right",
                        hidden_by_default=True,
                        visibility_group="ctf_micrograph",
                        width=4096,
                        height=3072,
                    ),
                ],
            }],
            "ctf_artifact_micrograph_toggle_available": True,
            "arguments": [],
            "logs": [],
            "auto_refresh": False,
        }

        rendered = render_to_string("nice_classic/batchview.html", context)

        self.assertNotIn("data-storage-key=", rendered)
        self.assertIn('id="batch_ctf_micrograph_toggle"', rendered)
        self.assertIn('aria-label="show micrographs"', rendered)
        self.assertIn("data-ctf-micrograph-preview", rendered)
        self.assertIn(
            'class="hidden relative flex-shrink-0 rounded-md overflow-hidden bg-black border border-streamdivider"',
            rendered,
        )
        self.assertIn("data-artifact-preview-group", rendered)
        self.assertEqual(rendered.count("data-artifact-preview-kind="), 2)
        self.assertNotIn("data-artifact-secondary-preview", rendered)
        self.assertIn('data-artifact-preview-kind="image"', rendered)
        self.assertIn('data-artifact-preview-kind="micrograph"', rendered)
        self.assertEqual(rendered.count("data-artifact-dimensions"), 2)
        self.assertIn("800 × 600 px", rendered)
        self.assertIn("4096 × 3072 px", rendered)
        self.assertGreater(rendered.index("movie_thumb.jpg"), rendered.index("movie_ctf_estimate_diag.jpg"))
        self.assertIn("setBatchCtfMicrographsVisible(false)", batch_view)
        self.assertIn('visible ? "hide micrographs" : "show micrographs"', batch_view)
        self.assertNotIn("nice_batch_ctf_micrograph_view_", batch_view)

    def test_batch_detail_motion_artifacts_toggle_composite_halves(self):
        batch_view = self._read_template("nice_classic/batchview.html")
        thumbnail_path = "/project/workspace/1_motion_correct/movie_thumb.jpg"
        context = {
            "jobid": 7,
            "artifact_counts": [{"extension": "JPG", "count": 1}],
            "artifact_images": [{
                "name": "movie_thumb.jpg",
                "path": thumbnail_path,
                "previews": [
                    self._preview(
                        thumbnail_path,
                        "movie_thumb.jpg",
                        "power_spectrum",
                        "left",
                        visibility_group="motion",
                        width=4096,
                        height=3072,
                    ),
                    self._preview(
                        thumbnail_path,
                        "movie_thumb.jpg",
                        "micrograph",
                        "right",
                        hidden_by_default=True,
                        visibility_group="motion",
                        width=4096,
                        height=3072,
                    ),
                ],
            }],
            "motion_artifact_toggle_available": True,
            "arguments": [],
            "logs": [],
            "auto_refresh": False,
        }

        rendered = render_to_string("nice_classic/batchview.html", context)

        self.assertIn('id="batch_motion_artifact_toggle"', rendered)
        self.assertIn('aria-label="show micrographs"', rendered)
        self.assertNotIn("data-storage-key=", rendered)
        self.assertEqual(rendered.count("data-artifact-preview-kind="), 2)
        self.assertIn('data-artifact-preview-kind="power_spectrum"', rendered)
        self.assertIn('data-artifact-preview-kind="micrograph"', rendered)
        self.assertIn('data-motion-artifact-preview="power_spectrum"', rendered)
        self.assertIn('data-motion-artifact-preview="micrograph"', rendered)
        self.assertEqual(rendered.count("data-artifact-dimensions"), 2)
        self.assertIn("4096 × 3072 px", rendered)
        self.assertIn('class="hidden relative', rendered)
        self.assertIn('setBatchMotionArtifactView("power_spectrum")', batch_view)
        self.assertIn('const showSideBySide = view === "side_by_side"', batch_view)
        self.assertIn('actionDescription = "show side by side"', batch_view)
        self.assertIn('actionDescription = "show power spectra"', batch_view)
        self.assertIn('currentView === "micrograph"', batch_view)
        self.assertIn('currentView === "side_by_side"', batch_view)
        self.assertNotIn("data-artifact-secondary-preview", rendered)

    def test_batch_pick_output_uses_artifact_grid_and_shared_center_overlay(self):
        batch_view = self._read_template("nice_classic/batchview.html")
        initial_pick = self._read_template("nice_stream/panelinitialpick.html")
        reference_pick = self._read_template("nice_stream/panelreferencepicking.html")
        slider = self._read_template("nice_stream/includes/_pick_micrograph_slider.html")
        slide = self._read_template("nice_stream/includes/_pick_micrograph_slide.html")
        batch_pick = self._read_template("nice_classic/includes/_batch_pick_micrograph_preview.html")
        slider_script = self._read_static("nice_lite/pick_micrograph_slider.js")
        shared_include = "nice_stream/includes/_pick_micrograph_slider.html"
        shared_slide = "nice_stream/includes/_pick_micrograph_slide.html"

        self.assertNotIn(shared_include, batch_view)
        self.assertIn(shared_include, initial_pick)
        self.assertIn(shared_include, reference_pick)
        self.assertIn(shared_slide, slider)
        self.assertIn(shared_slide, batch_pick)
        self.assertIn("nice_lite/pick_micrograph_slider.js", batch_view)
        self.assertIn("pick_micrograph_slider.js' %}?v=14", batch_view)
        self.assertIn("nice_lite/pick_micrograph_slider.js", initial_pick)
        self.assertIn("nice_lite/pick_micrograph_slider.js", reference_pick)
        self.assertIn("nice_stream/includes/_scroll_btn.html", slider)
        self.assertIn("data-pick-micrograph-overlay", batch_pick)
        self.assertIn("object-[100%]", slide)
        self.assertIn("context.arc(", slider_script)
        self.assertIn("context.strokeRect(", slider_script)
        self.assertIn('overlayMode === "boxes"', slider_script)
        self.assertIn('overlayMode === "circle"', slider_script)
        self.assertIn('button.dataset.pickOverlayModes || "points,circle,boxes"', slider_script)
        self.assertIn("width * overlayScale", slider_script)
        self.assertIn("const overlayScale = overlaySize / nativeSize", slider_script)
        self.assertIn("(overlaySize * scale) / 2", slider_script)
        self.assertIn("findNativeOverlaySize", slider_script)
        self.assertIn("Number(box.x) * scale", slider_script)
        self.assertIn("Number(box.y) * scale", slider_script)
        self.assertIn("--preprocess-contrast", slide)
        self.assertIn("--preprocess-brightness", slide)

        rendered = render_to_string("nice_classic/batchview.html", {
            "jobid": 7,
            "pick_micrographs": [{
                "path": "/project/workspace/1_motion/movie_thumb.jpg",
                "number": 9,
                "xdim": 4096,
                "ydim": 3072,
                "boxes": [{"x": 101, "y": 202, "width": 180, "height": 180}],
            }],
            "pick_box_overlay_available": True,
            "arguments": [],
            "logs": [],
            "artifact_counts": [],
            "artifact_images": [],
            "auto_refresh": False,
        })

        self.assertNotIn("pick locations", rendered)
        self.assertIn('id="batch_artifact_images" data-output-grid', rendered)
        self.assertIn('id="batch_output_tile_size" type="range"', rendered)
        self.assertIn('data-output-grid-target-id="batch_artifacts_panel"', rendered)
        self.assertIn("data-pick-micrograph-grid-preview", rendered)
        self.assertIn('id="batch_pick_overlay_toggle"', rendered)
        self.assertIn('data-pick-overlay-mode="points"', rendered)
        self.assertIn('data-pick-overlay-modes="points,circle,boxes"', rendered)
        self.assertIn('data-pick-overlay-name="pick"', rendered)
        self.assertIn("points overlay; show circle", rendered)
        self.assertIn("data-pick-overlay-toggle-label>points", rendered)
        self.assertIn('id="batch_pick_overlay_size" type="range"', rendered)
        self.assertLess(
            rendered.index('id="batch_pick_overlay_size"'),
            rendered.index('id="batch_pick_overlay_toggle"'),
        )
        self.assertIn('min="1" max="1000" step="1" value="100" disabled', rendered)
        self.assertIn('aria-label="overlay size in pixels"', rendered)
        self.assertIn("data-pick-overlay-size-control", rendered)
        self.assertIn('id="batch_pick_overlay_size_value" type="number"', rendered)
        self.assertIn("data-pick-overlay-size-number", rendered)
        self.assertIn('aria-label="enter overlay size in pixels"', rendered)
        self.assertIn('data-pick-overlay-size="0"', rendered)
        self.assertNotIn('type="color"', rendered)
        self.assertIn('id="batch_pick_overlay_color_control"', rendered)
        self.assertIn("data-pick-overlay-color-control", rendered)
        self.assertIn("data-pick-overlay-color-toggle", rendered)
        self.assertIn('aria-haspopup="dialog" aria-expanded="false"', rendered)
        self.assertIn('id="batch_pick_overlay_color_panel"', rendered)
        self.assertIn('role="dialog" aria-label="pick overlay color"', rendered)
        self.assertIn('id="batch_pick_overlay_color_disc"', rendered)
        self.assertIn("data-pick-overlay-color-disc-pointer", rendered)
        self.assertIn('width="192" height="192" tabindex="0"', rendered)
        self.assertIn("data-pick-overlay-color-presets", rendered)
        self.assertEqual(rendered.count("data-pick-overlay-color-preset="), 6)
        self.assertIn('data-pick-overlay-color-preset="#ff3b30"', rendered)
        self.assertIn('data-pick-overlay-color-preset="#ffd60a"', rendered)
        self.assertIn('data-pick-overlay-color-preset="#64dd17"', rendered)
        self.assertIn('data-pick-overlay-color-preset="#00e5ff"', rendered)
        self.assertIn('data-pick-overlay-color-preset="#60a5fa"', rendered)
        self.assertIn('data-pick-overlay-color-preset="#ff2dff"', rendered)
        red_preset = rendered.split(
            'data-pick-overlay-color-preset="#ff3b30"', 1,
        )[1].split("</button>", 1)[0]
        blue_preset = rendered.split(
            'data-pick-overlay-color-preset="#60a5fa"', 1,
        )[1].split("</button>", 1)[0]
        self.assertIn('aria-pressed="true"', red_preset)
        self.assertIn('aria-pressed="false"', blue_preset)
        self.assertIn('id="batch_pick_overlay_color_hex" type="text"', rendered)
        self.assertIn("data-pick-overlay-color-hex", rendered)
        self.assertIn('data-pick-overlay-color-channel="brightness"', rendered)
        self.assertNotIn("data-pick-overlay-color-toggle-value", rendered)
        self.assertNotIn('data-pick-overlay-color-channel="hue"', rendered)
        self.assertNotIn('data-pick-overlay-color-channel="saturation"', rendered)
        self.assertIn('current color #ff3b30', rendered)
        self.assertIn('data-pick-overlay-color-swatch aria-hidden="true"', rendered)
        self.assertIn('style="background-color:#ff3b30"', rendered)
        self.assertIn('value="#ff3b30" maxlength="7"', rendered)
        self.assertIn('data-pick-overlay-color="#ff3b30"', rendered)
        self.assertIn('data-pick-overlay-color-output="brightness"', rendered)
        self.assertIn('value="100"', rendered)
        self.assertIn('sizeInput.addEventListener("input"', slider_script)
        self.assertIn('sizeNumber.addEventListener("input"', slider_script)
        self.assertIn('sizeNumber.addEventListener("change"', slider_script)
        self.assertIn("sizeNumber.value = String(pixels)", slider_script)
        self.assertIn("sizeNumber.max = sizeInput.max", slider_script)
        self.assertIn("setSize(nativeSize || sizeInput.value)", slider_script)
        self.assertIn('sizeControl.classList.toggle("hidden", !adjustableSize)', slider_script)
        self.assertIn("if (sizeInput) sizeInput.disabled = !adjustableSize", slider_script)
        self.assertIn("if (sizeNumber) sizeNumber.disabled = !adjustableSize", slider_script)
        self.assertIn("normalizeHexColor", slider_script)
        self.assertIn('const defaultPickOverlayColor = "#ff3b30"', slider_script)
        self.assertIn("hexToHsv", slider_script)
        self.assertIn("hsvToRgb", slider_script)
        self.assertIn("hsvToHex", slider_script)
        self.assertIn("drawColorDisc", slider_script)
        self.assertIn("context.createImageData", slider_script)
        self.assertIn("context.putImageData", slider_script)
        self.assertIn("setColorFromDiscCoordinates", slider_script)
        self.assertIn("setColorPanelOpen", slider_script)
        self.assertIn("for (const preset of colorPresets)", slider_script)
        self.assertIn("preset.dataset.pickOverlayColorPreset", slider_script)
        self.assertIn('preset.addEventListener("click"', slider_script)
        self.assertIn('colorDisc?.addEventListener("pointerdown"', slider_script)
        self.assertIn('colorDisc?.addEventListener("pointermove"', slider_script)
        self.assertIn('colorDisc?.addEventListener("keydown"', slider_script)
        self.assertIn('event.key === "ArrowLeft"', slider_script)
        self.assertIn('event.key === "ArrowUp"', slider_script)
        self.assertIn('event.key === "Escape"', slider_script)
        self.assertIn('event.key === "Enter"', slider_script)
        self.assertIn('document.addEventListener("click"', slider_script)
        self.assertIn('input.addEventListener("input"', slider_script)
        self.assertIn("container.dataset.pickOverlayColor = color", slider_script)
        self.assertIn("overlayContainer?.dataset.pickOverlayColor", slider_script)
        self.assertNotIn("colorInput", slider_script)
        self.assertNotIn('overlayMode === "points"', slider_script)
        self.assertIn("pick_micrograph_slider.js?v=14", rendered)
        self.assertIn('data-pick-display-target-id="batch_artifact_images"', rendered)
        self.assertIn('id="batch_pick_display_min" type="range"', rendered)
        self.assertIn('min="0" max="254" step="1" value="0"', rendered)
        self.assertIn('id="batch_pick_display_max" type="range"', rendered)
        self.assertIn('min="1" max="255" step="1" value="255"', rendered)
        self.assertIn('id="batch_pick_display_brightness" type="range"', rendered)
        self.assertIn('min="-0.5" max="0.5" step="0.001" value="0"', rendered)
        self.assertIn('data-pick-display-output="brightness"', rendered)
        self.assertIn(">0.000</output>", rendered)
        self.assertIn('id="batch_pick_display_contrast" type="range"', rendered)
        self.assertIn('min="0.5" max="0.999" step="0.001" value="0.5"', rendered)
        self.assertIn('data-pick-display-output="contrast"', rendered)
        self.assertIn(">0.500</output>", rendered)
        self.assertIn("data-pick-display-reset", rendered)
        self.assertNotIn("are synchronized views", rendered)
        self.assertIn('const defaults = {min: 0, max: 255, brightness: 0, contrast: 0.5}', slider_script)
        self.assertIn("const applyLevels =", slider_script)
        self.assertIn("const applyBrightnessContrast =", slider_script)
        self.assertIn('name === "brightness" || name === "contrast"', slider_script)
        self.assertIn("minimum + maximum - (displayMinimum + displayMaximum)", slider_script)
        self.assertIn("displayRange * (1 - contrast)", slider_script)
        self.assertIn('target.style.setProperty("--pick-level-brightness"', slider_script)
        self.assertIn('target.style.setProperty("--pick-level-contrast"', slider_script)
        self.assertNotIn("--pick-display-brightness", slider_script)
        self.assertNotIn("--pick-display-contrast", slider_script)
        self.assertIn("pickDisplayImageSelector", slider_script)
        self.assertIn("target.querySelectorAll(imageSelector)", slider_script)
        self.assertIn("image.style.filter = displayFilter", slider_script)
        self.assertIn('querySelectorAll("[data-pick-display-controls]")', slider_script)
        self.assertNotIn("percent of picked box size", slider_script)
        self.assertNotIn("data-pick-micrograph-slider", rendered)
        self.assertIn('data-xdim="4096"', rendered)
        self.assertIn('data-ydim="3072"', rendered)
        self.assertIn("data-artifact-dimensions", rendered)
        self.assertIn("4096 × 3072 px", rendered)
        self.assertIn("movie_thumb.jpg", rendered)
        self.assertIn("picked micrograph 9", rendered)

    def test_extract_output_uses_lazy_paginated_particle_thumbnails(self):
        batch_view = self._read_template("nice_classic/batchview.html")
        context = {
            "jobid": 7,
            "particle_stack_page": {
                "stacks": [{"name": "particles.mrcs", "count": 481}],
                "particles": [
                    {
                        "number": 42,
                        "stack_name": "particles.mrcs",
                        "stack_index": 42,
                        "width": 128,
                        "height": 128,
                    },
                ],
                "total": 481,
                "page": 2,
                "pages": 12,
                "page_numbers": [1, 2, 3, 4, "…", 12],
                "ellipsis": "…",
                "has_previous": True,
                "previous_page": 1,
                "has_next": True,
                "next_page": 3,
                "first_particle": 41,
                "last_particle": 80,
            },
            "arguments": [],
            "logs": [],
            "artifact_counts": [{"extension": "MRCS", "count": 1}],
            "artifact_images": [],
            "auto_refresh": False,
        }

        rendered = render_to_string("nice_classic/batchview.html", context)

        self.assertIn("extracted particles", rendered)
        self.assertIn("41–80 of 481", rendered)
        self.assertIn("thumbnails are generated on demand", rendered)
        self.assertIn("/batchparticle/7/particles.mrcs/42", rendered)
        self.assertIn('loading="lazy" decoding="async"', rendered)
        self.assertIn('data-artifact-preview-kind="particle"', rendered)
        self.assertIn('id="batch_particle_thumbnails" data-output-grid', rendered)
        self.assertIn('id="batch_output_tile_size" type="range"', rendered)
        self.assertIn('data-output-grid-target-id="batch_artifacts_panel"', rendered)
        self.assertIn("data-output-grid-tile", rendered)
        self.assertIn("data-artifact-dimensions", rendered)
        self.assertIn("128 × 128 px", rendered)
        self.assertIn('style="height:var(--artifact-tile-height);width:var(--artifact-tile-height)"', rendered)
        self.assertIn("hover:opacity-90 transition-opacity", rendered)
        self.assertNotIn("<figcaption", rendered)
        self.assertNotIn(">#42<", rendered)
        self.assertIn('aria-label="extracted particle pages"', rendered)
        self.assertIn("?particle_page=1#batch_particle_gallery", rendered)
        self.assertIn("?particle_page=3#batch_particle_gallery", rendered)
        self.assertIn("?particle_page=12#batch_particle_gallery", rendered)
        self.assertIn('aria-label="last particle page"', rendered)
        self.assertIn("particle_stack_page.ellipsis", batch_view)
        self.assertIn("…", rendered)
        self.assertIn('aria-current="page"', rendered)
        self.assertNotIn("mrc2jpeg", batch_view)

    def test_import_movie_output_uses_lazy_paginated_thumbnail_tiles(self):
        batch_view = self._read_template("nice_classic/batchview.html")
        context = {
            "jobid": 7,
            "import_movie_page": {
                "movies": [{
                    "number": 41,
                    "name": "movie_041.mrcs",
                    "frames": 40,
                    "width": 4096,
                    "height": 4096,
                    "preview_available": True,
                    "token": "signed-token",
                }],
                "total": 481,
                "page": 2,
                "pages": 13,
                "page_numbers": [1, 2, 3, 4, "…", 13],
                "ellipsis": "…",
                "has_previous": True,
                "previous_page": 1,
                "has_next": True,
                "next_page": 3,
                "first_movie": 41,
                "last_movie": 80,
            },
            "arguments": [],
            "logs": [],
            "artifact_counts": [],
            "artifact_images": [],
            "auto_refresh": False,
        }

        rendered = render_to_string("nice_classic/batchview.html", context)

        self.assertIn("imported movies", rendered)
        self.assertIn("41–80 of 481", rendered)
        self.assertIn("thumbnails are cached as WebP on demand", rendered)
        self.assertIn("/batchmovie/7/signed-token", rendered)
        self.assertIn('loading="lazy" decoding="async"', rendered)
        self.assertIn('data-artifact-preview-kind="movie"', rendered)
        self.assertIn('id="batch_movie_thumbnails" data-output-grid', rendered)
        self.assertIn('id="batch_output_tile_size" type="range"', rendered)
        self.assertIn('data-output-grid-target-id="batch_artifacts_panel"', rendered)
        self.assertIn("data-output-grid-tile", rendered)
        self.assertIn("data-artifact-dimensions", rendered)
        self.assertIn("4096 × 4096 px", rendered)
        self.assertIn('style="height:var(--artifact-tile-height);width:var(--artifact-tile-height)"', rendered)
        self.assertNotIn("<figcaption", rendered)
        self.assertIn('aria-label="imported movie pages"', rendered)
        self.assertIn("?movie_page=1#batch_movie_gallery", rendered)
        self.assertIn("?movie_page=3#batch_movie_gallery", rendered)
        self.assertIn("?movie_page=13#batch_movie_gallery", rendered)
        self.assertIn('aria-label="last movie page"', rendered)
        self.assertIn("import_movie_page.ellipsis", batch_view)
        self.assertIn("…", rendered)

    def test_import_movie_output_shows_eer_dimensions_without_thumbnail(self):
        rendered = render_to_string(
            "nice_classic/batchview.html",
            {
                "jobid": 7,
                "import_movie_page": {
                    "movies": [{
                        "number": 1,
                        "name": "movie.eer",
                        "width": 4096,
                        "height": 4096,
                        "preview_available": False,
                        "token": "",
                    }],
                    "total": 1,
                    "page": 1,
                    "pages": 1,
                    "page_numbers": [1],
                    "ellipsis": "…",
                    "has_previous": False,
                    "has_next": False,
                    "first_movie": 1,
                    "last_movie": 1,
                },
                "arguments": [],
                "logs": [],
                "artifact_counts": [],
                "artifact_images": [],
                "auto_refresh": False,
            },
        )

        self.assertIn("movie.eer — 4096 × 4096 px", rendered)
        self.assertIn("data-artifact-dimensions", rendered)
        self.assertIn("thumbnail unavailable for this format", rendered)
        self.assertNotIn("/batchmovie/7/", rendered)

    def test_output_dimensions_use_page_local_header_switch(self):
        batch_view = self._read_template("nice_classic/batchview.html")
        rendered = render_to_string(
            "nice_classic/batchview.html",
            {
                "jobid": 7,
                "output_dimensions_available": True,
                "import_movie_page": {
                    "movies": [{
                        "number": 1,
                        "name": "movie.mrcs",
                        "width": 4096,
                        "height": 4096,
                        "preview_available": False,
                        "token": "",
                    }],
                    "total": 1,
                    "page": 1,
                    "pages": 1,
                    "page_numbers": [1],
                    "ellipsis": "…",
                    "has_previous": False,
                    "has_next": False,
                    "first_movie": 1,
                    "last_movie": 1,
                },
                "arguments": [],
                "logs": [],
                "artifact_counts": [],
                "artifact_images": [],
                "auto_refresh": False,
            },
        )

        self.assertIn('id="batch_output_controls"', rendered)
        self.assertIn('id="batch_output_tile_size" type="range"', rendered)
        self.assertIn('id="batch_output_dimensions_toggle"', rendered)
        self.assertLess(
            rendered.index('id="batch_output_tile_size"'),
            rendered.index('id="batch_output_dimensions_toggle"'),
        )
        self.assertIn('type="checkbox" role="switch"', rendered)
        self.assertIn('aria-label="hide output dimensions"', rendered)
        self.assertIn('class="peer sr-only" checked', rendered)
        self.assertIn('data-collapse-companion-id="batch_output_controls"', rendered)
        self.assertIn('querySelectorAll("[data-artifact-dimensions]")', batch_view)
        self.assertIn('dimensions.classList.toggle("hidden", !visible)', batch_view)
        self.assertIn("setBatchOutputDimensionsVisible(true)", batch_view)
        self.assertNotIn("nice_batch_output_dimensions", batch_view)

        rendered_without_dimensions = render_to_string(
            "nice_classic/batchview.html",
            {
                "jobid": 7,
                "output_dimensions_available": False,
                "arguments": [],
                "logs": [],
                "artifact_counts": [],
                "artifact_images": [],
                "auto_refresh": False,
            },
        )
        self.assertNotIn('id="batch_output_dimensions_toggle"', rendered_without_dimensions)

    def test_batch_detail_panels_collapse_and_logs_are_last(self):
        batch_view = self._read_template("nice_classic/batchview.html")

        for panel_key in ("overview", "arguments", "result", "artifacts", "logs"):
            self.assertIn(f'data-collapse-key="{panel_key}"', batch_view)
        self.assertNotIn('data-collapse-key="log-{{ forloop.counter }}"', batch_view)
        self.assertGreater(
            batch_view.index('data-collapse-key="logs"'),
            batch_view.index('data-collapse-key="artifacts"'),
        )
        self.assertEqual(batch_view.count('type="button" data-collapse-toggle'), 5)
        self.assertIn("data-collapse-icon", batch_view)
        self.assertIn("rotate-90 transition-transform", batch_view)
        self.assertIn('aria-expanded="true"', batch_view)
        self.assertIn("nice_batch_collapsed_panels_{{ jobid }}", batch_view)
        self.assertIn('classList.toggle("hidden", !expanded)', batch_view)
        self.assertIn('classList.toggle("rotate-90", expanded)', batch_view)

    def test_batch_detail_logs_support_tabs_and_split_view(self):
        batch_view = self._read_template("nice_classic/batchview.html")

        self.assertIn('role="tablist" aria-label="batch logs"', batch_view)
        self.assertIn('type="button" role="tab"', batch_view)
        self.assertIn('role="tabpanel"', batch_view)
        self.assertIn("data-log-tab", batch_view)
        self.assertIn("data-log-panel", batch_view)
        self.assertIn('id="batch_log_layout_toggle"', batch_view)
        self.assertIn("data-log-layout-label", batch_view)
        self.assertIn("nice_batch_log_tab_{{ jobid }}", batch_view)
        self.assertIn("nice_batch_log_layout_{{ jobid }}", batch_view)
        self.assertIn('window.matchMedia("(min-width: 64rem)")', batch_view)
        self.assertIn('batchLogWideScreenQuery.matches ? "split" : "tabs"', batch_view)
        self.assertIn('storedLayout === "split" || storedLayout === "tabs"', batch_view)
        self.assertIn("batchLogLayoutUserSelected", batch_view)
        self.assertIn('batchLogLayout === "split" ? "tabs" : "split"', batch_view)
        self.assertIn('event.key === "ArrowRight"', batch_view)
        self.assertIn('event.key === "ArrowLeft"', batch_view)
        self.assertIn('["grid", "grid-cols-1", "lg:grid-cols-3", "items-start", "gap-3"]', batch_view)
