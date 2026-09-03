(() => {
    const parseBoxes = (raw) => {
        try {
            return JSON.parse((raw || "[]").replaceAll("'", '"'));
        } catch (error) {
            return [];
        }
    };

    const drawPickOverlay = (preview, slide) => {
        const canvas = preview.querySelector("[data-pick-micrograph-overlay]");
        if (!canvas || !slide) return;

        const xdim = Number(slide.dataset.xdim || 0);
        const ydim = Number(slide.dataset.ydim || 0);
        const context = canvas.getContext("2d");
        context.clearRect(0, 0, canvas.width, canvas.height);
        if (!xdim || !ydim) return;

        const scale = xdim > ydim
            ? canvas.width / xdim
            : canvas.height / ydim;
        const xoffset = (canvas.width - (scale * xdim)) / 2;
        const yoffset = (canvas.height - (scale * ydim)) / 2;
        const markerColor = getComputedStyle(document.body)
            .getPropertyValue("--color-streamaction").trim() || "#60a5fa";
        const overlayMode = preview.closest("[data-pick-overlay-mode]")
            ?.dataset.pickOverlayMode || "points";

        context.strokeStyle = markerColor;
        for (const box of parseBoxes(slide.dataset.boxes)) {
            if (!box || !Number.isFinite(Number(box.x)) || !Number.isFinite(Number(box.y))) {
                continue;
            }
            const width = Number(box.width);
            const height = Number(box.height);
            if (
                overlayMode === "boxes"
                && Number.isFinite(width)
                && width > 0
                && Number.isFinite(height)
                && height > 0
            ) {
                context.strokeRect(
                    ((Number(box.x) - (width / 2)) * scale) + xoffset,
                    ((Number(box.y) - (height / 2)) * scale) + yoffset,
                    width * scale,
                    height * scale,
                );
                continue;
            }
            context.beginPath();
            context.arc(
                (Number(box.x) * scale) + xoffset,
                (Number(box.y) * scale) + yoffset,
                1,
                0,
                2 * Math.PI,
            );
            context.stroke();
        }
    };

    const sizePickOverlay = (preview) => {
        const canvas = preview.querySelector("[data-pick-micrograph-overlay]");
        const width = Math.round(preview.clientWidth);
        const height = Math.round(preview.clientHeight);
        if (!canvas || !width || !height) return;
        if (canvas.width !== width) canvas.width = width;
        if (canvas.height !== height) canvas.height = height;
    };

    const initializePreview = (preview) => {
        const slider = preview.querySelector("[data-pick-micrograph-slider]");
        const slides = [...preview.querySelectorAll("[data-pick-micrograph-slide]")];
        if (!slides.length) return;

        preview._activePickMicrographSlide = slides[0];
        sizePickOverlay(preview);
        drawPickOverlay(preview, slides[0]);
        if (typeof ResizeObserver === "function") {
            const resizeObserver = new ResizeObserver(() => {
                sizePickOverlay(preview);
                drawPickOverlay(preview, preview._activePickMicrographSlide);
            });
            resizeObserver.observe(preview);
        }
        if (!slider || typeof IntersectionObserver !== "function") return;

        const observer = new IntersectionObserver((entries) => {
            for (const entry of entries) {
                if (entry.intersectionRatio > 0.5) {
                    preview._activePickMicrographSlide = entry.target;
                    drawPickOverlay(preview, entry.target);
                }
            }
        }, {root: slider, threshold: 0.5});
        for (const slide of slides) observer.observe(slide);
    };

    const redrawPickPreviews = (container) => {
        container.querySelectorAll("[data-pick-micrograph-preview]").forEach((preview) => {
            const slide = preview._activePickMicrographSlide
                || preview.querySelector("[data-pick-micrograph-slide]");
            drawPickOverlay(preview, slide);
        });
    };

    const initializeOverlayToggle = (button) => {
        const container = document.getElementById(button.getAttribute("aria-controls"));
        if (!container) return;

        const setMode = (mode) => {
            const showBoxes = mode === "boxes";
            container.dataset.pickOverlayMode = showBoxes ? "boxes" : "points";
            button.setAttribute("aria-pressed", showBoxes ? "true" : "false");
            const actionDescription = showBoxes ? "show points" : "show boxes";
            button.setAttribute("aria-label", actionDescription);
            button.title = actionDescription;
            const label = button.querySelector("[data-pick-overlay-toggle-label]");
            if (label) label.textContent = actionDescription;
            redrawPickPreviews(container);
        };

        button.addEventListener("click", () => {
            const nextMode = container.dataset.pickOverlayMode === "boxes"
                ? "points"
                : "boxes";
            setMode(nextMode);
        });
        setMode("points");
    };

    window.scrollPickMicrographs = (button, event, direction) => {
        if (event) event.preventDefault();
        const carousel = button?.closest("[data-pick-micrograph-carousel]");
        const slider = carousel?.querySelector("[data-pick-micrograph-slider]");
        if (!slider) return;

        const forwards = Number(direction) > 0;
        const atStart = slider.scrollLeft <= 0;
        const atEnd = slider.scrollLeft + slider.clientWidth >= slider.scrollWidth - 1;
        if (forwards) {
            slider.scrollLeft = atEnd ? 0 : slider.scrollLeft + slider.clientWidth;
        } else {
            slider.scrollLeft = atStart
                ? slider.scrollWidth - slider.clientWidth
                : slider.scrollLeft - slider.clientWidth;
        }
    };

    const initializeAll = () => {
        document.querySelectorAll("[data-pick-micrograph-preview]").forEach(initializePreview);
        document.querySelectorAll("[data-pick-overlay-toggle]").forEach(initializeOverlayToggle);
    };
    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", initializeAll, {once: true});
    } else {
        initializeAll();
    }
})();
