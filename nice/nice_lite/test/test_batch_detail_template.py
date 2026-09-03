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
        self.assertIn('id="batch_artifact_images" class="flex flex-wrap items-start gap-2"', batch_view)
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
                    ),
                    self._preview(
                        micrograph_path,
                        "movie_thumb.jpg",
                        "micrograph",
                        "right",
                        hidden_by_default=True,
                        visibility_group="ctf_micrograph",
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
                    ),
                    self._preview(
                        thumbnail_path,
                        "movie_thumb.jpg",
                        "micrograph",
                        "right",
                        hidden_by_default=True,
                        visibility_group="motion",
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
        self.assertIn("nice_lite/pick_micrograph_slider.js", initial_pick)
        self.assertIn("nice_lite/pick_micrograph_slider.js", reference_pick)
        self.assertIn("nice_stream/includes/_scroll_btn.html", slider)
        self.assertIn("data-pick-micrograph-overlay", batch_pick)
        self.assertIn("object-[100%]", slide)
        self.assertIn("context.arc(", slider_script)
        self.assertIn("context.strokeRect(", slider_script)
        self.assertIn('overlayMode === "boxes"', slider_script)
        self.assertIn("Number(box.x) * scale", slider_script)
        self.assertIn("Number(box.y) * scale", slider_script)

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

        self.assertIn("pick locations", rendered)
        self.assertIn('id="batch_artifact_images" class="flex flex-wrap items-start gap-2"', rendered)
        self.assertIn("data-pick-micrograph-grid-preview", rendered)
        self.assertIn('id="batch_pick_overlay_toggle"', rendered)
        self.assertIn('data-pick-overlay-mode="points"', rendered)
        self.assertIn("show boxes", rendered)
        self.assertNotIn("data-pick-micrograph-slider", rendered)
        self.assertIn('data-xdim="4096"', rendered)
        self.assertIn('data-ydim="3072"', rendered)
        self.assertIn("movie_thumb.jpg", rendered)
        self.assertIn("picked micrograph 9", rendered)

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
