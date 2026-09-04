# global imports
import hashlib
import logging
import math
import os
import shutil
import signal
import tempfile
import time
from collections import Counter

import psutil
from PIL import Image

# django imports
from django.core.paginator import Paginator
from django.db import transaction
from django.utils import timezone

# local imports
from ..helpers import directory_exists, ensure_directory
from ..models import JobModel, WorkspaceModel
from .simple import SIMPLEBatch, SIMPLEProjFile, SIMPLEProject
from .job import Job
from .mrc import read_mrc_stack_info, render_mrc_particle_png
from .movie import render_movie_webp
from .workspace import Workspace


logger = logging.getLogger(__name__)


class BatchJob(Job):
    """Classic (non-stream) SIMPLE job attached to a workspace."""

    TERMINAL_STATUSES = frozenset(("finished", "failed", "stopped"))
    DELETABLE_STATUSES = TERMINAL_STATUSES | frozenset(("queued",))
    LOG_FILES = (
        ("stdout.log", "standard output"),
        ("stderr.log", "standard error"),
        ("nice_status.log", "NICE status callbacks"),
    )
    ARTIFACT_EXTENSIONS = frozenset((
        ".simple", ".mrc", ".star", ".jpg", ".jpeg", ".png", ".pdf", ".txt",
    ))
    IMAGE_EXTENSIONS = frozenset((".jpg", ".jpeg", ".png"))
    CTF_DIAGNOSTIC_SUFFIX = "_ctf_estimate_diag"
    MOTION_THUMBNAIL_SUFFIX = "_thumb.jpg"
    PICK_INTEGRATED_SUFFIX = "_intg"
    PARTICLE_STACK_PROGRAMS = frozenset(("extract", "reextract"))
    IMPORT_MOVIE_EXTENSIONS = frozenset((
        ".mrc", ".mrcs", ".tif", ".tiff", ".eer",
    ))
    MOVIE_THUMBNAIL_CACHE_DIR = ".nice_movie_thumbnails"
    MOVIE_THUMBNAIL_CACHE_VERSION = 2

    def __init__(self, pckg=None, id=None, request=None):
        super().__init__(id=None)
        self.disp = 0
        self.prnt = 0
        self.name = ""
        self.dirc = ""
        self.cdat = ""
        self.prog = ""
        self.args = {}
        self.wspc = None
        self.pckg = pckg
        self.status = "unknown"
        self.source = None

        if id is not None:
            self.id = id
            self.load()
        elif request is not None:
            self.setIDFromRequest(request)
            if self.id > 0:
                self.load()

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def setIDFromRequest(self, request):
        self.set_id_from_request(request, key="jobid")

    def load(self):
        """Populate fields from DB. Resets to empty state if not found."""
        self.jobmodel = JobModel.objects.filter(id=self.id).first()
        metadata = self.jobmodel.master_stats if self.jobmodel is not None else None
        if not isinstance(metadata, dict) or metadata.get("job_type") != "batch":
            self.jobmodel = None
            self.id = 0
            self.absdir = None
            return

        self.name = self.jobmodel.name
        self.desc = self.jobmodel.desc
        self.dirc = self.jobmodel.dirc
        self.cdat = self.jobmodel.cdat
        self.args = self.jobmodel.args
        self.wspc = self.jobmodel.dset
        self.disp = self.jobmodel.disp
        self.prog = metadata.get("program", "")
        self.pckg = metadata.get("package", "")
        self.prnt = metadata.get("parent", 0)
        self.status = self.jobmodel.status
        self.source = metadata.get("source")
        self.absdir = self.get_absdir()

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def get_absdir(self):
        if self.wspc is None:
            return None
        return os.path.join(self.wspc.proj.dirc, self.wspc.dirc, self.dirc)

    def getAbsDir(self):
        """Backward-compatible alias."""
        return self.get_absdir()

    def get_jobmodel(self):
        return self.jobmodel

    def get_safe_job_dir(self):
        """Return a validated, non-symlink job directory inside its workspace."""
        if self.jobmodel is None:
            return None

        project_root = os.path.realpath(self.jobmodel.dset.proj.dirc)
        workspace_raw = os.path.abspath(os.path.join(project_root, self.jobmodel.dset.dirc))
        workspace_root = os.path.realpath(workspace_raw)
        try:
            workspace_is_safe = os.path.commonpath((project_root, workspace_root)) == project_root
        except ValueError:
            workspace_is_safe = False
        if not workspace_is_safe or not os.path.isdir(workspace_root):
            return None

        job_dirc = self.jobmodel.dirc
        if (
            not isinstance(job_dirc, str)
            or job_dirc in ("", ".", "..")
            or job_dirc != os.path.basename(job_dirc)
        ):
            return None

        job_raw = os.path.abspath(os.path.join(workspace_root, job_dirc))
        job_root = os.path.realpath(job_raw)
        try:
            job_is_safe = (
                job_root != workspace_root
                and os.path.commonpath((workspace_root, job_root)) == workspace_root
            )
        except ValueError:
            job_is_safe = False
        if not job_is_safe or os.path.islink(job_raw) or not os.path.isdir(job_root):
            return None
        return job_root

    def _safe_job_file(self, filename, containment_root=None):
        """Return one fixed root-level job file after containment validation."""
        job_dir = self.get_safe_job_dir()
        if job_dir is None or not isinstance(filename, str) or filename != os.path.basename(filename):
            return None

        raw_path = os.path.abspath(os.path.join(job_dir, filename))
        resolved_path = os.path.realpath(raw_path)
        containment_root = os.path.realpath(containment_root or job_dir)
        try:
            path_is_safe = os.path.commonpath((containment_root, resolved_path)) == containment_root
        except ValueError:
            path_is_safe = False
        if not path_is_safe or not os.path.isfile(raw_path):
            return None
        return resolved_path

    def get_result_project_path(self):
        """Resolve this job's output project without trusting a request path."""
        job_dir = self.get_safe_job_dir()
        if job_dir is None or self.jobmodel is None:
            return None

        workspace_root = os.path.realpath(os.path.dirname(job_dir))
        workspace_project = self._safe_job_file("workspace.simple", workspace_root)
        if workspace_project is not None:
            return workspace_project

        if self.prog == "new_project":
            projname = str(self.args.get("projname", "")).strip()
            if projname and projname == os.path.basename(projname):
                project_filename = projname if projname.endswith(".simple") else f"{projname}.simple"
                named_project = self._safe_job_file(project_filename, workspace_root)
                if named_project is not None:
                    return named_project

        try:
            with os.scandir(job_dir) as directory_entries:
                filenames = sorted(
                    entry.name
                    for entry in directory_entries
                    if entry.name.endswith(".simple") and entry.is_file(follow_symlinks=True)
                )
        except OSError:
            return None

        candidates = []
        for filename in filenames:
            project_path = self._safe_job_file(filename, workspace_root)
            if project_path is not None and project_path not in candidates:
                candidates.append(project_path)
        return candidates[0] if len(candidates) == 1 else None

    def get_log_tails(self, max_bytes=131072):
        """Return bounded tails for the fixed batch log filenames."""
        job_dir = self.get_safe_job_dir()
        if job_dir is None:
            return []
        if not isinstance(max_bytes, int) or isinstance(max_bytes, bool) or max_bytes <= 0:
            max_bytes = 131072
        max_bytes = min(max_bytes, 1048576)

        logs = []
        for filename, label in self.LOG_FILES:
            path = self._safe_job_file(filename, job_dir)
            entry = {
                "name": filename,
                "label": label,
                "exists": path is not None,
                "size": 0,
                "truncated": False,
                "text": "",
            }
            if path is not None:
                try:
                    size = os.path.getsize(path)
                    with open(path, "rb") as log_file:
                        start = max(0, size - max_bytes)
                        log_file.seek(start)
                        log_text = log_file.read(max_bytes).decode("utf-8", errors="replace")
                    entry.update({
                        "size": size,
                        "truncated": start > 0,
                        "text": log_text,
                    })
                except OSError:
                    entry["exists"] = False
            logs.append(entry)
        return logs

    def _get_source_batch_job(self):
        """Return this job's earlier, same-workspace batch source."""
        if self.jobmodel is None:
            return None

        metadata = self.jobmodel.master_stats
        source = metadata.get("source") if isinstance(metadata, dict) else None
        if not isinstance(source, dict):
            return None
        source_job_id = source.get("batch_job_id")
        if (
            source.get("type") != "batch_job"
            or not isinstance(source_job_id, int)
            or isinstance(source_job_id, bool)
            or source_job_id <= 0
        ):
            return None

        source_model = JobModel.objects.filter(
            id=source_job_id,
            dset_id=self.jobmodel.dset_id,
        ).first()
        source_metadata = source_model.master_stats if source_model is not None else None
        if (
            source_model is None
            or source_model.disp >= self.jobmodel.disp
            or not isinstance(source_metadata, dict)
            or source_metadata.get("job_type") != "batch"
        ):
            return None
        return BatchJob(id=source_job_id)

    def _get_source_motion_job(self):
        """Return an earlier motion job from this job's persisted source chain."""
        if self.jobmodel is None or self.prog not in ("ctf_estimate", "pick"):
            return None

        source_job = self
        visited_job_ids = {self.id}
        while source_job is not None:
            source_job = source_job._get_source_batch_job()
            if source_job is None or source_job.id in visited_job_ids:
                return None
            visited_job_ids.add(source_job.id)
            if source_job.prog == "motion_correct":
                return source_job
        return None

    def _get_motion_thumbnail(self, source_motion_job, diagnostic_name):
        """Resolve the source thumbnail matching one CTF diagnostic filename."""
        if source_motion_job is None:
            return None
        diagnostic_stem, extension = os.path.splitext(diagnostic_name)
        if (
            extension.lower() not in self.IMAGE_EXTENSIONS
            or not diagnostic_stem.endswith(self.CTF_DIAGNOSTIC_SUFFIX)
        ):
            return None

        micrograph_stem = diagnostic_stem[:-len(self.CTF_DIAGNOSTIC_SUFFIX)]
        if not micrograph_stem:
            return None
        source_dir = source_motion_job.get_safe_job_dir()
        if source_dir is None:
            return None
        return source_motion_job._safe_job_file(
            f"{micrograph_stem}{self.MOTION_THUMBNAIL_SUFFIX}",
            source_dir,
        )

    def _get_motion_micrograph_dimensions(self, thumbnail_name):
        """Return dimensions of the integrated MRC represented by a thumbnail."""
        if (
            self.prog != "motion_correct"
            or not isinstance(thumbnail_name, str)
            or not thumbnail_name.lower().endswith(self.MOTION_THUMBNAIL_SUFFIX)
        ):
            return None
        motion_dir = self.get_safe_job_dir()
        if motion_dir is None:
            return None
        micrograph_stem = thumbnail_name[:-len(self.MOTION_THUMBNAIL_SUFFIX)]
        integrated_path = self._safe_job_file(
            f"{micrograph_stem}{self.PICK_INTEGRATED_SUFFIX}.mrc",
            motion_dir,
        )
        return (
            self._read_mrc_dimensions(integrated_path)
            if integrated_path is not None
            else None
        )

    @staticmethod
    def _read_mrc_dimensions(path):
        """Read validated x/y dimensions through the shared MRC helper."""
        stack_info = read_mrc_stack_info(path)
        if stack_info is None:
            return None
        return stack_info.width, stack_info.height

    @staticmethod
    def _read_raster_dimensions(path):
        """Return JPEG/PNG dimensions without loading image pixels."""
        try:
            with Image.open(path) as image:
                width, height = image.size
        except (
            EOFError,
            Image.DecompressionBombError,
            OSError,
            OverflowError,
            TypeError,
            ValueError,
        ):
            return None
        if (
            not isinstance(width, int)
            or isinstance(width, bool)
            or not isinstance(height, int)
            or isinstance(height, bool)
            or width <= 0
            or height <= 0
        ):
            return None
        return width, height

    @staticmethod
    def _read_box_centers(path, max_coordinates):
        """Read bounded EMAN box records and convert top-left corners to centers."""
        centers = []
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as box_file:
                for line_number, line in enumerate(box_file):
                    if line_number >= max_coordinates:
                        break
                    fields = line.split()
                    if len(fields) < 4:
                        continue
                    try:
                        xcoord, ycoord, width, height = map(float, fields[:4])
                    except ValueError:
                        continue
                    if (
                        not all(math.isfinite(value) for value in (xcoord, ycoord, width, height))
                        or width <= 0
                        or height <= 0
                    ):
                        continue
                    centers.append({
                        "x": xcoord + (width / 2.0),
                        "y": ycoord + (height / 2.0),
                        "width": width,
                        "height": height,
                    })
        except OSError:
            return []
        return centers

    def get_pick_micrograph_previews(self, max_previews=20, max_coordinates=1500):
        """Build stream-shaped picker previews from owned batch artifacts."""
        if self.jobmodel is None or self.prog != "pick":
            return []
        if not isinstance(max_previews, int) or isinstance(max_previews, bool) or max_previews <= 0:
            max_previews = 20
        if (
            not isinstance(max_coordinates, int)
            or isinstance(max_coordinates, bool)
            or max_coordinates <= 0
        ):
            max_coordinates = 1500
        max_previews = min(max_previews, 24)
        max_coordinates = min(max_coordinates, 5000)

        pick_dir = self.get_safe_job_dir()
        source_motion_job = self._get_source_motion_job()
        source_dir = (
            source_motion_job.get_safe_job_dir()
            if source_motion_job is not None
            else None
        )
        if pick_dir is None or source_dir is None:
            return []

        try:
            with os.scandir(pick_dir) as directory_entries:
                box_names = sorted(
                    entry.name
                    for entry in directory_entries
                    if entry.name.lower().endswith(".box")
                    and entry.is_file(follow_symlinks=True)
                )
        except OSError:
            return []

        selected_box_names = box_names[-max_previews:]
        first_number = len(box_names) - len(selected_box_names) + 1
        previews = []
        for number, box_name in enumerate(selected_box_names, start=first_number):
            box_path = self._safe_job_file(box_name, pick_dir)
            box_stem = os.path.splitext(box_name)[0]
            if box_path is None or not box_stem.endswith(self.PICK_INTEGRATED_SUFFIX):
                continue

            micrograph_stem = box_stem[:-len(self.PICK_INTEGRATED_SUFFIX)]
            thumbnail_path = source_motion_job._safe_job_file(
                f"{micrograph_stem}{self.MOTION_THUMBNAIL_SUFFIX}",
                source_dir,
            )
            integrated_path = source_motion_job._safe_job_file(
                f"{box_stem}.mrc",
                source_dir,
            )
            if thumbnail_path is None or integrated_path is None:
                continue
            dimensions = self._read_mrc_dimensions(integrated_path)
            if dimensions is None:
                continue

            previews.append({
                "path": thumbnail_path,
                "number": number,
                "xdim": dimensions[0],
                "ydim": dimensions[1],
                "boxes": self._read_box_centers(box_path, max_coordinates),
            })
        return previews

    @staticmethod
    def _artifact_preview(
        path,
        name,
        kind="image",
        crop="full",
        hidden_by_default=False,
        visibility_group=None,
        dimensions=None,
    ):
        """Return one normalized image view for the batch artifact gallery."""
        if kind == "power_spectrum":
            suffix = "power spectrum"
        elif kind == "micrograph":
            suffix = "micrograph"
        else:
            suffix = ""
        preview = {
            "path": path,
            "name": name,
            "kind": kind,
            "crop": crop,
            "title": f"{name} — {suffix}" if suffix else name,
            "alt": f"{name} {suffix}" if suffix else name,
        }
        if hidden_by_default:
            preview["hidden_by_default"] = True
        if visibility_group is not None:
            preview["visibility_group"] = visibility_group
        if dimensions is not None:
            preview["width"], preview["height"] = dimensions
        return preview

    def get_artifact_summary(self, max_previews=12):
        """Summarize recognized root-level result files and safe image previews."""
        job_dir = self.get_safe_job_dir()
        if job_dir is None:
            return {"counts": [], "images": []}
        if not isinstance(max_previews, int) or isinstance(max_previews, bool) or max_previews <= 0:
            max_previews = 12
        max_previews = min(max_previews, 24)

        try:
            with os.scandir(job_dir) as directory_entries:
                entries = sorted(directory_entries, key=lambda entry: entry.name)
        except OSError:
            return {"counts": [], "images": []}

        source_motion_job = self._get_source_motion_job()

        counts = Counter()
        images = []
        for entry in entries:
            extension = os.path.splitext(entry.name)[1].lower()
            if extension not in self.ARTIFACT_EXTENSIONS:
                continue
            path = self._safe_job_file(entry.name, job_dir)
            if path is None:
                continue
            counts[extension] += 1
            if extension in self.IMAGE_EXTENSIONS and len(images) < max_previews:
                raster_dimensions = self._read_raster_dimensions(path)
                image = {
                    "name": entry.name,
                    "path": path,
                    "previews": [self._artifact_preview(
                        path,
                        entry.name,
                        dimensions=raster_dimensions,
                    )],
                }
                if (
                    self.prog == "motion_correct"
                    and entry.name.lower().endswith(self.MOTION_THUMBNAIL_SUFFIX)
                ):
                    micrograph_dimensions = (
                        self._get_motion_micrograph_dimensions(entry.name)
                        or raster_dimensions
                    )
                    image["previews"] = [
                        self._artifact_preview(
                            path,
                            entry.name,
                            kind="power_spectrum",
                            crop="left",
                            visibility_group="motion",
                            dimensions=micrograph_dimensions,
                        ),
                        self._artifact_preview(
                            path,
                            entry.name,
                            kind="micrograph",
                            crop="right",
                            hidden_by_default=True,
                            visibility_group="motion",
                            dimensions=micrograph_dimensions,
                        ),
                    ]
                micrograph_path = self._get_motion_thumbnail(source_motion_job, entry.name)
                if micrograph_path is not None:
                    micrograph_name = os.path.basename(micrograph_path)
                    micrograph_dimensions = (
                        source_motion_job._get_motion_micrograph_dimensions(
                            micrograph_name
                        )
                        or source_motion_job._read_raster_dimensions(micrograph_path)
                    )
                    image["previews"].append(self._artifact_preview(
                        micrograph_path,
                        micrograph_name,
                        kind="micrograph",
                        crop="right",
                        hidden_by_default=True,
                        visibility_group="ctf_micrograph",
                        dimensions=micrograph_dimensions,
                    ))
                images.append(image)

        return {
            "counts": [
                {"extension": extension.removeprefix(".").upper(), "count": count}
                for extension, count in sorted(counts.items())
            ],
            "images": images,
        }

    def get_particle_stack_page(self, page=1, page_size=40):
        """Return one page of addressable images from owned extract stacks.

        Only MRC headers are read here. Pixel data is read later by the image
        endpoint for the thumbnails that the browser actually requests.
        """
        empty_page = {
            "stacks": [],
            "particles": [],
            "total": 0,
            "page": 1,
            "pages": 0,
            "page_numbers": [],
            "ellipsis": Paginator.ELLIPSIS,
            "has_previous": False,
            "has_next": False,
            "first_particle": 0,
            "last_particle": 0,
        }
        if self.prog not in self.PARTICLE_STACK_PROGRAMS:
            return empty_page
        if not isinstance(page, int) or isinstance(page, bool) or page < 1:
            page = 1
        if not isinstance(page_size, int) or isinstance(page_size, bool) or page_size < 1:
            page_size = 40
        page_size = min(page_size, 100)

        job_dir = self.get_safe_job_dir()
        if job_dir is None:
            return empty_page
        try:
            with os.scandir(job_dir) as directory_entries:
                stack_names = sorted(
                    entry.name
                    for entry in directory_entries
                    if entry.name.lower().endswith(".mrcs")
                    and entry.is_file(follow_symlinks=True)
                )
        except OSError:
            return empty_page

        stacks = []
        total = 0
        for stack_name in stack_names:
            stack_path = self._safe_job_file(stack_name, job_dir)
            stack_info = read_mrc_stack_info(stack_path) if stack_path is not None else None
            if stack_info is None:
                continue
            stacks.append({
                "name": stack_name,
                "count": stack_info.count,
                "width": stack_info.width,
                "height": stack_info.height,
                "first_particle": total + 1,
            })
            total += stack_info.count

        if total == 0:
            return empty_page

        paginator = Paginator(range(total), page_size)
        page_obj = paginator.get_page(page)
        first_offset = page_obj.start_index() - 1
        last_offset = page_obj.end_index()
        particles = []
        stack_offset = 0
        for stack in stacks:
            stack_last_offset = stack_offset + stack["count"]
            selected_start = max(first_offset, stack_offset)
            selected_end = min(last_offset, stack_last_offset)
            for global_offset in range(selected_start, selected_end):
                particles.append({
                    "number": global_offset + 1,
                    "stack_name": stack["name"],
                    "stack_index": global_offset - stack_offset + 1,
                    "width": stack["width"],
                    "height": stack["height"],
                })
            stack_offset = stack_last_offset
            if stack_offset >= last_offset:
                break

        return {
            "stacks": stacks,
            "particles": particles,
            "total": total,
            "page": page_obj.number,
            "pages": paginator.num_pages,
            "page_numbers": list(paginator.get_elided_page_range(
                page_obj.number,
                on_each_side=2,
                on_ends=1,
            )),
            "ellipsis": Paginator.ELLIPSIS,
            "has_previous": page_obj.has_previous(),
            "previous_page": (
                page_obj.previous_page_number() if page_obj.has_previous() else None
            ),
            "has_next": page_obj.has_next(),
            "next_page": page_obj.next_page_number() if page_obj.has_next() else None,
            "first_particle": page_obj.start_index(),
            "last_particle": page_obj.end_index(),
        }

    def get_particle_thumbnail(self, stack_name, particle_index, max_size=160):
        """Return one owned extract particle as in-memory PNG bytes."""
        if (
            self.prog not in self.PARTICLE_STACK_PROGRAMS
            or not isinstance(stack_name, str)
            or stack_name != os.path.basename(stack_name)
            or not stack_name.lower().endswith(".mrcs")
        ):
            return None
        job_dir = self.get_safe_job_dir()
        stack_path = self._safe_job_file(stack_name, job_dir) if job_dir is not None else None
        if stack_path is None:
            return None
        return render_mrc_particle_png(stack_path, particle_index, max_size=max_size)

    def get_import_movie_thumbnail(self, movie_path, max_size=160):
        """Return one imported movie thumbnail from an on-demand WebP cache."""
        if (
            self.prog != "import_movies"
            or not isinstance(movie_path, str)
            or not os.path.isabs(movie_path)
            or not os.path.isfile(movie_path)
        ):
            return None

        cache_path = self._import_movie_thumbnail_cache_path(movie_path, max_size)
        if cache_path is None:
            return None
        try:
            with open(cache_path, "rb") as cache_file:
                return cache_file.read()
        except FileNotFoundError:
            pass
        except OSError:
            return None

        thumbnail = render_movie_webp(movie_path, max_size=max_size)
        if thumbnail is None:
            return None

        temporary_path = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="wb",
                dir=os.path.dirname(cache_path),
                prefix=".movie-thumbnail-",
                suffix=".tmp",
                delete=False,
            ) as temporary_file:
                temporary_path = temporary_file.name
                temporary_file.write(thumbnail)
            os.replace(temporary_path, cache_path)
            temporary_path = None
        except OSError:
            return thumbnail
        finally:
            if temporary_path is not None:
                try:
                    os.unlink(temporary_path)
                except OSError:
                    pass
        return thumbnail

    def get_import_movie_paths(self):
        """Return import candidates from this job's submitted movie source."""
        if self.prog != "import_movies" or not isinstance(self.args, dict):
            return []

        directory = str(self.args.get("dir_movies", "") or "").strip()
        file_table = str(self.args.get("filetab", "") or "").strip()
        if bool(directory) == bool(file_table):
            return []

        job_dir = self.get_safe_job_dir()
        if job_dir is None:
            return []
        if directory:
            if not os.path.isabs(directory):
                directory = os.path.abspath(os.path.join(job_dir, directory))
            try:
                with os.scandir(directory) as entries:
                    return sorted(
                        os.path.abspath(entry.path)
                        for entry in entries
                        if entry.is_file(follow_symlinks=True)
                        and os.path.splitext(entry.name)[1].lower()
                        in self.IMPORT_MOVIE_EXTENSIONS
                    )
            except OSError:
                return []

        if not os.path.isabs(file_table):
            file_table = os.path.abspath(os.path.join(job_dir, file_table))
        try:
            with open(file_table, encoding="utf-8") as source_file:
                source_lines = source_file.readlines()
        except OSError:
            return []

        movie_paths = []
        for source_line in source_lines:
            movie_path = source_line.strip()
            if not movie_path or movie_path.startswith("#"):
                continue
            if not os.path.isabs(movie_path):
                movie_path = os.path.abspath(os.path.join(job_dir, movie_path))
            if (
                os.path.splitext(movie_path)[1].lower()
                in self.IMPORT_MOVIE_EXTENSIONS
                and os.path.isfile(movie_path)
            ):
                movie_paths.append(movie_path)
        return movie_paths

    def _import_movie_thumbnail_cache_path(self, movie_path, max_size):
        """Return a safe cache path keyed by source identity and dimensions."""
        if (
            not isinstance(max_size, int)
            or isinstance(max_size, bool)
            or max_size < 1
        ):
            return None
        job_dir = self.get_safe_job_dir()
        if job_dir is None:
            return None
        cache_dir = os.path.join(job_dir, self.MOVIE_THUMBNAIL_CACHE_DIR)
        if os.path.islink(cache_dir) or not ensure_directory(cache_dir):
            return None
        try:
            if os.path.commonpath((job_dir, os.path.realpath(cache_dir))) != job_dir:
                return None
            source_stat = os.stat(movie_path)
        except (OSError, ValueError):
            return None

        cache_identity = "\0".join((
            str(self.MOVIE_THUMBNAIL_CACHE_VERSION),
            os.path.realpath(movie_path),
            str(source_stat.st_size),
            str(source_stat.st_mtime_ns),
            str(max_size),
        ))
        cache_name = hashlib.sha256(cache_identity.encode("utf-8")).hexdigest()
        return os.path.join(cache_dir, cache_name + ".webp")

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _create_dir(self, parent_dir):
        if self.dirc == "":
            logger.error("_create_dir: empty dir name")
            return False
        if not directory_exists(parent_dir):
            logger.error("_create_dir: parent directory missing")
            return False

        new_dir_path = os.path.join(parent_dir, self.dirc)
        if directory_exists(new_dir_path):
            logger.error("_create_dir: destination already exists")
            return False

        return ensure_directory(new_dir_path)

    def createDir(self, parent_dir):
        """Backward-compatible alias."""
        return self._create_dir(parent_dir)

    def createLink(self, source, destination):
        if not os.path.exists(source):
            logger.error("createLink: source missing")
            return False
        if not os.path.isfile(source):
            logger.error("createLink: source is not a file")
            return False
        try:
            os.symlink(source, destination)
        except OSError:
            logger.error("createLink: symlink creation failed")
            return False
        return True

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    @staticmethod
    def _metadata(pckg, prog, parent=0, source=None):
        metadata = {
            "job_type": "batch",
            "package": pckg,
            "program": prog,
            "parent": parent,
        }
        if isinstance(source, dict):
            metadata["source"] = dict(source)
        return metadata

    def new(self, workspace, pckg, prog, args, parent_proj=None, source=None, display_name=None):
        """Create and launch a SIMPLE or SINGLE batch job."""
        if workspace is None or pckg not in ("simple", "single") or not prog:
            logger.error("new: invalid batch job configuration")
            return False
        if not isinstance(args, dict):
            logger.error("new: args is not a dictionary")
            return False

        workspacemodel = workspace.get_workspacemodel()
        workspace_dir = workspace.get_absdir()
        if workspacemodel is None or not directory_exists(workspace_dir):
            logger.error("new: workspace is unavailable")
            return False

        explicit_parent_proj = parent_proj is not None
        if explicit_parent_proj:
            try:
                workspace_root = os.path.realpath(workspace_dir)
                parent_proj = os.path.realpath(parent_proj)
                source_in_workspace = os.path.commonpath((workspace_root, parent_proj)) == workspace_root
            except (TypeError, ValueError):
                source_in_workspace = False
            if not source_in_workspace or not os.path.isfile(parent_proj):
                logger.error("new: batch project source is unavailable")
                return False
        else:
            parent_proj = os.path.join(workspace_dir, "workspace.simple")
            if not os.path.isfile(parent_proj):
                if not SIMPLEProject(workspace_dir).create():
                    logger.error("new: failed to initialize workspace.simple")
                    return False

        self.pckg = pckg
        self.prog = prog
        if isinstance(display_name, str) and display_name.strip():
            self.name = display_name.strip()
        else:
            self.name = prog.replace("_", " ")
        self.args = dict(args)

        # Reserve the display counter under a row lock so two near-simultaneous
        # Start clicks cannot choose the same job directory.
        with transaction.atomic():
            locked_workspace = WorkspaceModel.objects.select_for_update().get(pk=workspacemodel.pk)
            self.disp = locked_workspace.jcnt + 1
            self.dirc = f"{self.disp}_{prog}"
            jobmodel = JobModel(
                dset=locked_workspace,
                cdat=timezone.now(),
                disp=self.disp,
                args=self.args,
                name=self.name,
                dirc=self.dirc,
                status="queued",
                master_status="queued",
                master_stats=self._metadata(pckg, prog, source=source),
                master_heartbeat=0,
            )
            jobmodel.save()
            locked_workspace.jcnt = self.disp
            locked_workspace.save()

        self.id = jobmodel.id
        self.jobmodel = jobmodel
        self.wspc = locked_workspace
        self.status = "queued"
        self.absdir = self.get_absdir()

        if not self._create_dir(workspace_dir):
            self.status = "failed"
            jobmodel.status = "failed"
            jobmodel.master_status = "failed"
            jobmodel.save()
            return False

        simple = SIMPLEBatch(pckg=pckg)
        launch_options = {"parent_proj": parent_proj} if explicit_parent_proj else {}
        if simple.start(self.args, self.absdir, workspace_dir, prog, self.id, **launch_options):
            return True

        self.status = "failed"
        jobmodel.status = "failed"
        jobmodel.master_status = "failed"
        jobmodel.save()
        return False

    def linkParticleSet(self, project, workspace, set_proj):
        self.args = {}

        workspacemodel = WorkspaceModel.objects.filter(id=workspace.id).first()
        if workspacemodel is None:
            logger.error("linkParticleSet: workspace not found")
            return False

        self.disp = workspacemodel.jcnt + 1
        self.name = "particle set"
        self.dirc = str(self.disp) + "_link_particle_set"
        workspace_dir = os.path.join(project.dirc, workspacemodel.dirc)
        if not self._create_dir(workspace_dir):
            return False

        link_dst = os.path.join(workspace_dir, self.dirc, "workspace.simple")
        if not self.createLink(set_proj, link_dst):
            return False

        jobmodel = JobModel(
            dset=workspacemodel,
            cdat=timezone.now(),
            disp=self.disp,
            args={},
            name=self.name,
            dirc=self.dirc,
            status="finished",
            master_status="finished",
            master_stats=self._metadata("simple_stream", "link_particle_set"),
        )
        jobmodel.save()
        workspacemodel.jcnt = self.disp
        workspacemodel.save()

        self.id = jobmodel.id
        self.jobmodel = jobmodel
        self.wspc = workspacemodel
        self.status = "finished"
        self.absdir = self.get_absdir()
        return True

    def linkParticleSetFinal(self, project, workspace, set_proj, set_desel):
        self.args = {}
        self.prnt = 0

        workspacemodel = WorkspaceModel.objects.filter(id=workspace.id).first()
        if workspacemodel is None:
            logger.error("linkParticleSetFinal: workspace not found")
            return False

        self.disp = workspacemodel.jcnt + 1
        self.pckg = "simple"
        self.prog = "selection"
        self.name = "particle set"
        self.dirc = str(self.disp) + "_" + self.prog
        workspace_dir = os.path.join(project.dirc, workspacemodel.dirc)
        if not self._create_dir(workspace_dir):
            return False

        self.args["deselfile"] = set_desel
        self.args["oritype"] = "cls2D"

        jobmodel = JobModel(
            dset=workspacemodel,
            cdat=timezone.now(),
            disp=self.disp,
            name=self.name,
            args=self.args,
            dirc=self.dirc,
            status="queued",
            master_status="queued",
            master_stats=self._metadata(self.pckg, self.prog, self.prnt),
        )
        jobmodel.save()
        workspacemodel.jcnt = self.disp
        workspacemodel.save()

        self.id = jobmodel.id
        self.jobmodel = jobmodel
        self.wspc = workspacemodel
        self.status = "queued"
        self.absdir = self.get_absdir()

        simple = SIMPLEBatch(pckg=self.pckg)
        return simple.start(
            self.args,
            os.path.join(workspace_dir, self.dirc),
            workspace_dir,
            self.prog,
            self.id,
            parent_proj=set_proj,
        )

    # ------------------------------------------------------------------
    # Updates / completion
    # ------------------------------------------------------------------

    def _local_process_is_running(self, pid):
        """Return True while the generated local job process still exists."""
        try:
            process = psutil.Process(pid)
            if not process.is_running() or process.status() == psutil.STATUS_ZOMBIE:
                return False
            process_dir = os.path.realpath(process.cwd())
        except (psutil.NoSuchProcess, psutil.ZombieProcess):
            return False
        except psutil.AccessDenied:
            return True
        return process_dir == os.path.realpath(self.get_absdir() or "")

    def queued_job_can_delete(self):
        """Return True for a stale local queued job whose process has exited."""
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is None or jobmodel.status != "queued":
            return False

        job_dir = self.get_absdir()
        if job_dir is None:
            return False

        try:
            with open(os.path.join(job_dir, "job.script"), encoding="utf-8") as script_file:
                if "# localtemplate" not in script_file.read(4096):
                    return False
            with open(os.path.join(job_dir, "nice.pid"), encoding="utf-8") as pid_file:
                pid = int(pid_file.read(32).strip())
        except (OSError, ValueError):
            return False

        return pid > 1 and not self._local_process_is_running(pid)

    def reconcile_local_completion(self):
        """Mark a locally dispatched job finished after a verified normal exit.

        This is a conservative fallback for local jobs whose final NICE callback
        did not reach Django. Scheduler jobs and abnormal local exits are left
        unchanged for their authoritative status path to handle.
        """
        if self.jobmodel is None or self.status not in ("queued", "running"):
            return False

        job_dir = self.get_absdir()
        if job_dir is None:
            return False

        try:
            with open(os.path.join(job_dir, "job.script"), encoding="utf-8") as script_file:
                if "# localtemplate" not in script_file.read(4096):
                    return False
            with open(os.path.join(job_dir, "nice.pid"), encoding="utf-8") as pid_file:
                pid = int(pid_file.read(32).strip())
        except (OSError, ValueError):
            return False

        if pid <= 1 or self._local_process_is_running(pid):
            return False

        try:
            stdout_path = os.path.join(job_dir, "stdout.log")
            stdout_size = os.path.getsize(stdout_path)
            with open(stdout_path, "rb") as stdout_file:
                stdout_file.seek(max(0, stdout_size - 65536))
                stdout_tail = stdout_file.read().decode("utf-8", errors="replace")
        except OSError:
            return False

        normal_stop = any(
            line.startswith("**** ") and line.endswith(" NORMAL STOP ****")
            for line in stdout_tail.splitlines()
        )
        if not normal_stop:
            return False

        with transaction.atomic():
            jobmodel = JobModel.objects.select_for_update().filter(id=self.id).first()
            if jobmodel is None or jobmodel.status not in ("queued", "running"):
                return False
            jobmodel.status = "finished"
            jobmodel.master_status = "finished"
            jobmodel.save(update_fields=("status", "master_status"))

        self.jobmodel = jobmodel
        self.status = "finished"
        return True

    def _get_local_process_group(self, pid):
        """Return an isolated process group owned by this job, or ``None``."""
        job_dir = self.get_absdir()
        if job_dir is None:
            logger.error("stop: job directory is unavailable")
            return None

        try:
            process = psutil.Process(pid)
            if not process.is_running() or process.status() == psutil.STATUS_ZOMBIE:
                raise psutil.NoSuchProcess(pid)
            process_dir = os.path.realpath(process.cwd())
            job_dir = os.path.realpath(job_dir)
            process_group = os.getpgid(pid)
        except (OSError, psutil.Error):
            logger.error("stop: job process is not available on this host")
            return None

        if process_dir != job_dir:
            logger.error("stop: pid does not belong to the job directory")
            return None
        if process_group != pid:
            logger.error("stop: job process group is not isolated")
            return None
        return process_group

    def stop(self):
        """Stop a running local batch job by terminating its process group."""
        with transaction.atomic():
            jobmodel = JobModel.objects.select_for_update().filter(id=self.id).first()
            if jobmodel is None or jobmodel.status != "running":
                logger.error("stop: batch job is not running")
                return False

            pid_path = os.path.join(self.get_absdir() or "", "nice.pid")
            try:
                with open(pid_path, encoding="utf-8") as pid_file:
                    pid = int(pid_file.read(32).strip())
            except (OSError, ValueError):
                logger.error("stop: nice.pid is missing or invalid")
                return False
            if pid <= 1:
                logger.error("stop: nice.pid contains an unsafe process id")
                return False

            process_group = self._get_local_process_group(pid)
            if process_group is None:
                return False
            try:
                os.killpg(process_group, signal.SIGTERM)
            except OSError:
                logger.error("stop: failed to terminate the job process group")
                return False

            self.status = "stopped"
            jobmodel.status = "stopped"
            jobmodel.master_status = "stopped"
            jobmodel.master_update = {}
            jobmodel.save(update_fields=("status", "master_status", "master_update"))
            return True

    def updateDescription(self, request):
        if "new_job_description" in request.POST:
            return self.set_description(request.POST["new_job_description"])
        return False

    def markComplete(self, project, workspace):
        del project, workspace
        self.status = "finished"
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return False
        jobmodel.status = self.status
        jobmodel.master_status = self.status
        jobmodel.save()
        return True

    def updateStats(self, stats_json, project, workspace):
        del project, workspace
        with transaction.atomic():
            jobmodel = JobModel.objects.select_for_update().filter(id=self.id).first()
            if jobmodel is None:
                return {}
            # Script-level and in-process callbacks can race at shutdown. Once
            # one of them records a terminal state, do not let a later callback
            # resurrect the job or replace the original outcome.
            if jobmodel.status in self.TERMINAL_STATUSES:
                return {}

            response = {}
            if "job_heartbeat" in stats_json:
                jobmodel.master_heartbeat = int(time.time())

            job_stats = stats_json.get("job")
            if isinstance(job_stats, dict) and "terminate" in job_stats:
                terminal_status = job_stats.get("status")
                if not isinstance(terminal_status, str) or terminal_status not in self.TERMINAL_STATUSES:
                    terminal_status = "finished"
                jobmodel.status = terminal_status
                jobmodel.master_status = terminal_status
            else:
                jobmodel.status = "running"
                jobmodel.master_status = "running"
                response = jobmodel.master_update or {}
                jobmodel.master_update = {}

            jobmodel.save()
            return response

    # ------------------------------------------------------------------
    # Projfile helpers
    # ------------------------------------------------------------------

    def getProjectStats(self):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return None
        projfile = os.path.join(jobmodel.dset.proj.dirc, jobmodel.dset.dirc, self.dirc, "workspace.simple")
        return SIMPLEProjFile(projfile).getGlobalStats()

    def getProjectFieldStats(self, oritype, fromp=None, top=None, sortkey=None, sortasc=None, hist=False, boxes=False, plotkey=None):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return None
        projfile = os.path.join(jobmodel.dset.proj.dirc, jobmodel.dset.dirc, self.dirc, "workspace.simple")
        return SIMPLEProjFile(projfile).getFieldStats(oritype, fromp, top, sortkey, sortasc, hist, boxes, plotkey)

    # ------------------------------------------------------------------
    # Deletion
    # ------------------------------------------------------------------

    def delete(self, project, workspace):
        """Permanently remove a batch job and reset its workspace job counter."""
        del project
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return False
        if jobmodel.status not in self.DELETABLE_STATUSES:
            logger.error("delete: batch job status is not deletable")
            return False
        if jobmodel.status == "queued" and not self.queued_job_can_delete():
            logger.error("delete: queued batch job is still active or unverifiable")
            return False

        workspace_dir = os.path.realpath(
            os.path.join(jobmodel.dset.proj.dirc, jobmodel.dset.dirc)
        )
        if not directory_exists(workspace_dir):
            return False

        job_dirc = jobmodel.dirc
        if (
            not isinstance(job_dirc, str)
            or not job_dirc.strip()
            or os.path.isabs(job_dirc)
        ):
            return False

        raw_job_path = os.path.abspath(os.path.join(workspace_dir, job_dirc))
        job_path = os.path.realpath(raw_job_path)
        try:
            path_is_safe = (
                job_path != workspace_dir
                and os.path.commonpath((workspace_dir, job_path)) == workspace_dir
            )
        except ValueError:
            path_is_safe = False
        if not path_is_safe or os.path.islink(raw_job_path):
            logger.error("delete: unsafe batch job directory")
            return False
        if not directory_exists(job_path):
            return False

        try:
            shutil.rmtree(job_path)
        except OSError:
            logger.error("delete: failed to permanently remove batch job directory")
            return False

        workspace_id = jobmodel.dset_id
        jobmodel.delete()
        if (
            getattr(workspace, "id", None) != workspace_id
            or not callable(getattr(workspace, "reconcile_job_counter", None))
        ):
            workspace = Workspace(workspace_id)
        workspace.reconcile_job_counter()
        return True
