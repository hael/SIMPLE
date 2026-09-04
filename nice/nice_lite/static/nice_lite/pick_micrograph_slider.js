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
        const overlayContainer = preview.closest("[data-pick-overlay-mode]");
        const overlayMode = overlayContainer?.dataset.pickOverlayMode || "points";
        const requestedOverlaySize = Number(
            overlayContainer?.dataset.pickOverlaySize || 0,
        );

        context.strokeStyle = markerColor;
        for (const box of parseBoxes(slide.dataset.boxes)) {
            if (!box || !Number.isFinite(Number(box.x)) || !Number.isFinite(Number(box.y))) {
                continue;
            }
            const width = Number(box.width);
            const height = Number(box.height);
            if (
                Number.isFinite(width)
                && width > 0
                && Number.isFinite(height)
                && height > 0
            ) {
                const nativeSize = Math.min(width, height);
                const overlaySize = Number.isFinite(requestedOverlaySize)
                    && requestedOverlaySize > 0
                    ? requestedOverlaySize
                    : nativeSize;
                const overlayScale = overlaySize / nativeSize;
                const scaledWidth = width * overlayScale;
                const scaledHeight = height * overlayScale;
                if (overlayMode === "boxes") {
                    context.strokeRect(
                        ((Number(box.x) - (scaledWidth / 2)) * scale) + xoffset,
                        ((Number(box.y) - (scaledHeight / 2)) * scale) + yoffset,
                        scaledWidth * scale,
                        scaledHeight * scale,
                    );
                    continue;
                }
                if (overlayMode === "circle") {
                    context.beginPath();
                    context.arc(
                        (Number(box.x) * scale) + xoffset,
                        (Number(box.y) * scale) + yoffset,
                        (overlaySize * scale) / 2,
                        0,
                        2 * Math.PI,
                    );
                    context.stroke();
                    continue;
                }
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

    const findNativeOverlaySize = (container) => {
        for (const slide of container.querySelectorAll("[data-pick-micrograph-slide]")) {
            for (const box of parseBoxes(slide.dataset.boxes)) {
                const width = Number(box?.width);
                const height = Number(box?.height);
                if (
                    Number.isFinite(width)
                    && width > 0
                    && Number.isFinite(height)
                    && height > 0
                ) {
                    return Math.round(Math.min(width, height));
                }
            }
        }
        return null;
    };

    const initializeOverlayToggle = (button) => {
        const container = document.getElementById(button.getAttribute("aria-controls"));
        if (!container) return;
        const modes = ["points", "circle", "boxes"];
        const sizeInput = document.getElementById(
            button.dataset.pickOverlaySizeInputId || "",
        );
        const sizeControl = sizeInput?.closest("[data-pick-overlay-size-control]");
        const sizeNumber = sizeControl?.querySelector("[data-pick-overlay-size-number]");

        const setSize = (rawPixels) => {
            if (!sizeInput) return;
            const minimum = Number(sizeInput.min) || 1;
            const maximum = Number(sizeInput.max) || 1000;
            const requestedPixels = Number(rawPixels);
            const pixels = Number.isFinite(requestedPixels)
                ? Math.round(Math.min(maximum, Math.max(minimum, requestedPixels)))
                : 100;
            sizeInput.value = String(pixels);
            sizeInput.setAttribute(
                "aria-valuetext",
                `${pixels} pixels`,
            );
            if (sizeNumber) {
                sizeNumber.value = String(pixels);
                sizeNumber.setAttribute("aria-valuetext", `${pixels} pixels`);
            }
            container.dataset.pickOverlaySize = String(pixels);
            redrawPickPreviews(container);
        };

        const setMode = (mode) => {
            const activeMode = modes.includes(mode) ? mode : "points";
            const nextMode = modes[(modes.indexOf(activeMode) + 1) % modes.length];
            container.dataset.pickOverlayMode = activeMode;
            button.dataset.pickOverlayMode = activeMode;
            button.setAttribute(
                "aria-label",
                `${activeMode} overlay; show ${nextMode}`,
            );
            button.title = `show ${nextMode}`;
            const label = button.querySelector("[data-pick-overlay-toggle-label]");
            if (label) label.textContent = activeMode;
            const adjustableSize = activeMode === "circle" || activeMode === "boxes";
            if (sizeControl) sizeControl.classList.toggle("hidden", !adjustableSize);
            if (sizeInput) sizeInput.disabled = !adjustableSize;
            if (sizeNumber) sizeNumber.disabled = !adjustableSize;
            redrawPickPreviews(container);
        };

        button.addEventListener("click", () => {
            const activeMode = modes.includes(button.dataset.pickOverlayMode)
                ? button.dataset.pickOverlayMode
                : "points";
            const nextMode = modes[(modes.indexOf(activeMode) + 1) % modes.length];
            setMode(nextMode);
        });
        if (sizeInput) {
            sizeInput.addEventListener("input", () => setSize(sizeInput.value));
            const nativeSize = findNativeOverlaySize(container);
            if (nativeSize) {
                sizeInput.max = String(Math.max(
                    Number(sizeInput.max) || 1000,
                    nativeSize * 3,
                ));
            }
            if (sizeNumber) {
                sizeNumber.min = sizeInput.min;
                sizeNumber.max = sizeInput.max;
                sizeNumber.addEventListener("input", () => {
                    if (sizeNumber.value !== "") setSize(sizeNumber.value);
                });
                sizeNumber.addEventListener("change", () => {
                    setSize(sizeNumber.value || sizeInput.value);
                });
            }
            setSize(nativeSize || sizeInput.value);
        }
        setMode("points");
    };

    const initializeDisplayControls = (controls) => {
        const target = document.getElementById(controls.dataset.pickDisplayTargetId || "");
        if (!target) return;

        const displayFilter = [
            "brightness(var(--pick-level-brightness,1))",
            "contrast(var(--pick-level-contrast,1))",
        ].join(" ");
        for (const image of target.querySelectorAll("[data-pick-micrograph-slide] > img")) {
            image.style.filter = displayFilter;
        }

        const inputs = Object.fromEntries(
            [...controls.querySelectorAll("[data-pick-display-control]")].map((input) => [
                input.dataset.pickDisplayControl,
                input,
            ]),
        );
        const outputs = Object.fromEntries(
            [...controls.querySelectorAll("[data-pick-display-output]")].map((output) => [
                output.dataset.pickDisplayOutput,
                output,
            ]),
        );
        if (!inputs.min || !inputs.max || !inputs.brightness || !inputs.contrast) return;

        const displayMinimum = 0;
        const displayMaximum = 255;
        const displayRange = displayMaximum - displayMinimum;
        const displayMidpoint = (displayMinimum + displayMaximum) / 2;
        const defaults = {min: 0, max: 255, brightness: 0, contrast: 0.5};
        const clamp = (value, minimum, maximum) => Math.min(
            maximum,
            Math.max(minimum, Number(value) || 0),
        );
        const setInputValue = (name, value) => {
            inputs[name].value = String(value);
        };
        const updateValue = (name, value, ariaValue = String(value)) => {
            if (outputs[name]) outputs[name].textContent = String(value);
            inputs[name].setAttribute("aria-valuetext", ariaValue);
        };

        const applyLevels = (rawMinimum, rawMaximum, changedName = "") => {
            let minimum = Math.round(clamp(rawMinimum, displayMinimum, displayMaximum - 1));
            let maximum = Math.round(clamp(rawMaximum, displayMinimum + 1, displayMaximum));
            if (minimum >= maximum) {
                if (changedName === "min") {
                    maximum = Math.min(displayMaximum, minimum + 1);
                } else {
                    minimum = Math.max(displayMinimum, maximum - 1);
                }
            }

            // EMAN2's inspector exposes brightness/contrast as an alternate
            // representation of the same min/max density window.
            const brightness = -0.5 * (
                minimum + maximum - (displayMinimum + displayMaximum)
            ) / displayRange;
            const contrast = 1 - (
                (minimum - maximum) / (2 * (displayMinimum - displayMaximum))
            );
            const formattedBrightness = brightness.toFixed(3);
            const formattedContrast = contrast.toFixed(3);

            setInputValue("min", minimum);
            setInputValue("max", maximum);
            setInputValue("brightness", formattedBrightness);
            setInputValue("contrast", formattedContrast);
            updateValue("min", minimum, `${minimum} of ${displayMaximum}`);
            updateValue("max", maximum, `${maximum} of ${displayMaximum}`);
            updateValue("brightness", formattedBrightness);
            updateValue("contrast", formattedContrast);

            const normalizedMinimum = minimum / displayMaximum;
            const normalizedMaximum = maximum / displayMaximum;
            const levelSum = normalizedMinimum + normalizedMaximum;
            const levelRange = normalizedMaximum - normalizedMinimum;
            target.style.setProperty("--pick-level-brightness", String(1 / levelSum));
            target.style.setProperty("--pick-level-contrast", String(levelSum / levelRange));
        };

        const applyBrightnessContrast = () => {
            const brightness = clamp(inputs.brightness.value, -0.5, 0.5);
            const contrast = clamp(inputs.contrast.value, 0.5, 0.999);
            const minimum = displayMidpoint
                - (displayRange * (1 - contrast))
                - (brightness * displayRange);
            const maximum = displayMidpoint
                + (displayRange * (1 - contrast))
                - (brightness * displayRange);
            applyLevels(minimum, maximum);
        };

        for (const [name, input] of Object.entries(inputs)) {
            input.addEventListener("input", () => {
                if (name === "brightness" || name === "contrast") {
                    applyBrightnessContrast();
                } else {
                    applyLevels(inputs.min.value, inputs.max.value, name);
                }
            });
        }
        controls.querySelector("[data-pick-display-reset]")?.addEventListener("click", () => {
            for (const [name, value] of Object.entries(defaults)) {
                setInputValue(name, value);
            }
            applyLevels(defaults.min, defaults.max);
        });
        applyLevels(defaults.min, defaults.max);
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
        document.querySelectorAll("[data-pick-display-controls]").forEach(initializeDisplayControls);
    };
    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", initializeAll, {once: true});
    } else {
        initializeAll();
    }
})();
