(() => {
    const defaultPickOverlayColor = "#ff3b30";

    const parseBoxes = (raw) => {
        try {
            return JSON.parse((raw || "[]").replaceAll("'", '"'));
        } catch (error) {
            return [];
        }
    };

    const normalizeHexColor = (rawColor) => {
        let color = String(rawColor || "").trim();
        if (!color.startsWith("#")) color = `#${color}`;
        return /^#[0-9a-f]{6}$/i.test(color) ? color.toLowerCase() : null;
    };

    const hexToHsv = (rawColor) => {
        const color = normalizeHexColor(rawColor) || defaultPickOverlayColor;
        const red = Number.parseInt(color.slice(1, 3), 16) / 255;
        const green = Number.parseInt(color.slice(3, 5), 16) / 255;
        const blue = Number.parseInt(color.slice(5, 7), 16) / 255;
        const maximum = Math.max(red, green, blue);
        const minimum = Math.min(red, green, blue);
        const delta = maximum - minimum;
        let hue = 0;
        if (delta) {
            if (maximum === red) hue = 60 * (((green - blue) / delta) % 6);
            else if (maximum === green) hue = 60 * (((blue - red) / delta) + 2);
            else hue = 60 * (((red - green) / delta) + 4);
        }
        if (hue < 0) hue += 360;
        return {
            hue: Math.round(hue),
            saturation: Math.round(maximum ? (delta / maximum) * 100 : 0),
            brightness: Math.round(maximum * 100),
        };
    };

    const hsvToRgb = (rawHue, rawSaturation, rawBrightness) => {
        const hue = ((Number(rawHue) % 360) + 360) % 360;
        const saturation = Math.min(100, Math.max(0, Number(rawSaturation))) / 100;
        const brightness = Math.min(100, Math.max(0, Number(rawBrightness))) / 100;
        const chroma = brightness * saturation;
        const intermediate = chroma * (1 - Math.abs(((hue / 60) % 2) - 1));
        let [red, green, blue] = [0, 0, 0];
        if (hue < 60) [red, green, blue] = [chroma, intermediate, 0];
        else if (hue < 120) [red, green, blue] = [intermediate, chroma, 0];
        else if (hue < 180) [red, green, blue] = [0, chroma, intermediate];
        else if (hue < 240) [red, green, blue] = [0, intermediate, chroma];
        else if (hue < 300) [red, green, blue] = [intermediate, 0, chroma];
        else [red, green, blue] = [chroma, 0, intermediate];
        const offset = brightness - chroma;
        return [red, green, blue].map((channel) => Math.round((channel + offset) * 255));
    };

    const hsvToHex = (rawHue, rawSaturation, rawBrightness) => {
        const channelToHex = (channel) => channel.toString(16).padStart(2, "0");
        return `#${hsvToRgb(rawHue, rawSaturation, rawBrightness)
            .map(channelToHex).join("")}`;
    };

    const overlayModeLabel = (mode) => mode === "boxes" ? "box" : mode;

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
        const defaultMarkerColor = getComputedStyle(document.body)
            .getPropertyValue("--color-streamaction").trim() || "#60a5fa";
        const overlayContainer = preview.closest("[data-pick-overlay-mode]");
        const overlayMode = overlayContainer?.dataset.pickOverlayMode || "points";
        const markerColor = overlayContainer?.dataset.pickOverlayColor
            || defaultMarkerColor;
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

    const sizeClassAverageOverlay = (preview) => {
        const canvas = preview.querySelector("[data-class-average-overlay]");
        const width = Math.round(preview.clientWidth);
        const height = Math.round(preview.clientHeight);
        if (!canvas || !width || !height) return;
        if (canvas.width !== width) canvas.width = width;
        if (canvas.height !== height) canvas.height = height;
    };

    const drawClassAverageOverlay = (preview) => {
        const canvas = preview.querySelector("[data-class-average-overlay]");
        if (!canvas) return;

        const sourceWidth = Number(canvas.dataset.sourceWidth || 0);
        const sourceHeight = Number(canvas.dataset.sourceHeight || 0);
        const context = canvas.getContext("2d");
        context.clearRect(0, 0, canvas.width, canvas.height);
        if (!sourceWidth || !sourceHeight) return;

        const overlayContainer = preview.closest("[data-pick-overlay-mode]");
        const overlayMode = overlayContainer?.dataset.pickOverlayMode || "circle";
        const requestedOverlaySize = Number(
            overlayContainer?.dataset.pickOverlaySize || 0,
        );
        const nativeSize = Math.min(sourceWidth, sourceHeight);
        const overlaySize = Number.isFinite(requestedOverlaySize)
            && requestedOverlaySize > 0
            ? requestedOverlaySize
            : nativeSize;
        const scale = Math.min(
            canvas.width / sourceWidth,
            canvas.height / sourceHeight,
        );
        const scaledOverlaySize = overlaySize * scale;
        const centerX = canvas.width / 2;
        const centerY = canvas.height / 2;

        context.strokeStyle = overlayContainer?.dataset.pickOverlayColor
            || defaultPickOverlayColor;
        if (overlayMode === "boxes") {
            context.strokeRect(
                centerX - (scaledOverlaySize / 2),
                centerY - (scaledOverlaySize / 2),
                scaledOverlaySize,
                scaledOverlaySize,
            );
        } else if (overlayMode === "circle") {
            context.beginPath();
            context.arc(centerX, centerY, scaledOverlaySize / 2, 0, 2 * Math.PI);
            context.stroke();
        }
    };

    let classAverageResizeObserver = null;

    const initializeClassAverageOverlay = (preview) => {
        sizeClassAverageOverlay(preview);
        drawClassAverageOverlay(preview);
        if (typeof ResizeObserver !== "function") return;
        if (!classAverageResizeObserver) {
            classAverageResizeObserver = new ResizeObserver((entries) => {
                for (const entry of entries) {
                    sizeClassAverageOverlay(entry.target);
                    drawClassAverageOverlay(entry.target);
                }
            });
        }
        classAverageResizeObserver.observe(preview);
    };

    const redrawOverlays = (container) => {
        container.querySelectorAll("[data-pick-micrograph-preview]").forEach((preview) => {
            const slide = preview._activePickMicrographSlide
                || preview.querySelector("[data-pick-micrograph-slide]");
            drawPickOverlay(preview, slide);
        });
        container.querySelectorAll("[data-class-average-overlay-preview]")
            .forEach(drawClassAverageOverlay);
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
        const classOverlay = container.querySelector("[data-class-average-overlay]");
        const sourceWidth = Number(classOverlay?.dataset.sourceWidth);
        const sourceHeight = Number(classOverlay?.dataset.sourceHeight);
        if (
            Number.isFinite(sourceWidth)
            && sourceWidth > 0
            && Number.isFinite(sourceHeight)
            && sourceHeight > 0
        ) {
            return Math.round(Math.min(sourceWidth, sourceHeight));
        }
        return null;
    };

    const initializeOverlayToggle = (button) => {
        const container = document.getElementById(button.getAttribute("aria-controls"));
        if (!container) return;
        const modes = (button.dataset.pickOverlayModes || "points,circle,boxes")
            .split(",")
            .map((mode) => mode.trim())
            .filter((mode) => ["points", "circle", "boxes"].includes(mode));
        const defaultMode = modes[0] || "points";
        const overlayName = button.dataset.pickOverlayName || "pick";
        const sizeInput = document.getElementById(
            button.dataset.pickOverlaySizeInputId || "",
        );
        const sizeControl = sizeInput?.closest("[data-pick-overlay-size-control]");
        const sizeNumber = sizeControl?.querySelector("[data-pick-overlay-size-number]");
        const colorControl = document.getElementById(
            button.dataset.pickOverlayColorControlId || "",
        );
        const colorToggle = colorControl?.querySelector("[data-pick-overlay-color-toggle]");
        const colorPanel = colorControl?.querySelector("[data-pick-overlay-color-panel]");
        const colorHex = colorControl?.querySelector("[data-pick-overlay-color-hex]");
        const colorDisc = colorControl?.querySelector("[data-pick-overlay-color-disc]");
        const colorDiscPointer = colorControl?.querySelector(
            "[data-pick-overlay-color-disc-pointer]",
        );
        const colorSwatches = colorControl?.querySelectorAll(
            "[data-pick-overlay-color-swatch]",
        ) || [];
        const colorPresets = colorControl?.querySelectorAll(
            "[data-pick-overlay-color-preset]",
        ) || [];
        const colorSliders = Object.fromEntries(
            [...(colorControl?.querySelectorAll("[data-pick-overlay-color-channel]") || [])]
                .map((input) => [input.dataset.pickOverlayColorChannel, input]),
        );
        const colorOutputs = Object.fromEntries(
            [...(colorControl?.querySelectorAll("[data-pick-overlay-color-output]") || [])]
                .map((output) => [output.dataset.pickOverlayColorOutput, output]),
        );
        let colorState = hexToHsv(
            container.dataset.pickOverlayColor || defaultPickOverlayColor,
        );
        let renderedDiscBrightness = null;
        let colorDiscDragging = false;

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
            redrawOverlays(container);
        };

        const drawColorDisc = () => {
            if (!colorDisc || renderedDiscBrightness === colorState.brightness) return;
            const context = colorDisc.getContext("2d");
            const {width, height} = colorDisc;
            const centerX = (width - 1) / 2;
            const centerY = (height - 1) / 2;
            const radius = Math.min(width, height) / 2;
            const image = context.createImageData(width, height);
            for (let y = 0; y < height; y += 1) {
                for (let x = 0; x < width; x += 1) {
                    const deltaX = x - centerX;
                    const deltaY = y - centerY;
                    const distance = Math.hypot(deltaX, deltaY);
                    const offset = ((y * width) + x) * 4;
                    if (distance > radius) {
                        image.data[offset + 3] = 0;
                        continue;
                    }
                    const hue = (Math.atan2(deltaY, deltaX) * 180 / Math.PI + 360) % 360;
                    const saturation = Math.min(100, (distance / radius) * 100);
                    const [red, green, blue] = hsvToRgb(
                        hue,
                        saturation,
                        colorState.brightness,
                    );
                    image.data[offset] = red;
                    image.data[offset + 1] = green;
                    image.data[offset + 2] = blue;
                    image.data[offset + 3] = Math.round(
                        Math.min(1, Math.max(0, radius - distance + 0.5)) * 255,
                    );
                }
            }
            context.putImageData(image, 0, 0);
            renderedDiscBrightness = colorState.brightness;
        };

        const updateColorDiscPointer = () => {
            if (!colorDisc || !colorDiscPointer) return;
            const width = colorDisc.clientWidth || colorDisc.width;
            const height = colorDisc.clientHeight || colorDisc.height;
            const radius = Math.min(width, height) / 2;
            const angle = colorState.hue * Math.PI / 180;
            const distance = radius * colorState.saturation / 100;
            colorDiscPointer.style.left = `${(width / 2) + (Math.cos(angle) * distance)}px`;
            colorDiscPointer.style.top = `${(height / 2) + (Math.sin(angle) * distance)}px`;
            colorDisc.setAttribute(
                "aria-valuetext",
                `hue ${colorState.hue} degrees, saturation ${colorState.saturation} percent`,
            );
        };

        const updateColorInterface = (color) => {
            for (const swatch of colorSwatches) swatch.style.backgroundColor = color;
            for (const preset of colorPresets) {
                preset.setAttribute(
                    "aria-pressed",
                    normalizeHexColor(preset.dataset.pickOverlayColorPreset) === color
                        ? "true"
                        : "false",
                );
            }
            if (colorToggle) {
                colorToggle.setAttribute(
                    "aria-label",
                    `open ${overlayName} overlay color picker, current color ${color}`,
                );
            }
            if (colorHex) {
                colorHex.value = color;
                colorHex.setAttribute("aria-invalid", "false");
            }
            for (const channel of ["hue", "saturation", "brightness"]) {
                const value = colorState[channel];
                if (colorSliders[channel]) {
                    colorSliders[channel].value = String(value);
                    colorSliders[channel].setAttribute(
                        "aria-valuetext",
                        channel === "hue" ? `${value} degrees` : `${value} percent`,
                    );
                }
                if (colorOutputs[channel]) {
                    colorOutputs[channel].textContent = channel === "hue"
                        ? `${value}°`
                        : `${value}%`;
                }
            }
            if (colorSliders.brightness) {
                const brightnessEnd = hsvToHex(
                    colorState.hue,
                    colorState.saturation,
                    100,
                );
                colorSliders.brightness.style.background =
                    `linear-gradient(to right,#000000,${brightnessEnd})`;
            }
            drawColorDisc();
            updateColorDiscPointer();
        };

        const commitColor = (color) => {
            container.dataset.pickOverlayColor = color;
            updateColorInterface(color);
            redrawOverlays(container);
        };

        const setColor = (rawColor) => {
            const color = normalizeHexColor(rawColor);
            if (!color) return false;
            colorState = hexToHsv(color);
            commitColor(color);
            return true;
        };

        const setColorFromDiscCoordinates = (clientX, clientY) => {
            if (!colorDisc) return;
            const bounds = colorDisc.getBoundingClientRect();
            const radius = Math.min(bounds.width, bounds.height) / 2;
            if (!radius) return;
            const deltaX = clientX - bounds.left - (bounds.width / 2);
            const deltaY = clientY - bounds.top - (bounds.height / 2);
            const distance = Math.min(radius, Math.hypot(deltaX, deltaY));
            if (distance > 0) {
                colorState.hue = Math.round(
                    (Math.atan2(deltaY, deltaX) * 180 / Math.PI + 360) % 360,
                ) % 360;
            }
            colorState.saturation = Math.round((distance / radius) * 100);
            commitColor(hsvToHex(
                colorState.hue,
                colorState.saturation,
                colorState.brightness,
            ));
        };

        const setColorPanelOpen = (open) => {
            if (!colorToggle || !colorPanel) return;
            colorPanel.classList.toggle("hidden", !open);
            colorToggle.setAttribute("aria-expanded", open ? "true" : "false");
            if (open) {
                drawColorDisc();
                updateColorDiscPointer();
                (colorDisc || colorHex)?.focus();
            }
        };

        const setMode = (mode) => {
            const activeMode = modes.includes(mode) ? mode : defaultMode;
            const nextMode = modes[(modes.indexOf(activeMode) + 1) % modes.length];
            const activeModeLabel = overlayModeLabel(activeMode);
            const nextModeLabel = overlayModeLabel(nextMode);
            container.dataset.pickOverlayMode = activeMode;
            button.dataset.pickOverlayMode = activeMode;
            button.setAttribute(
                "aria-label",
                `${activeModeLabel} overlay; show ${nextModeLabel}`,
            );
            button.title = `show ${nextModeLabel}`;
            const label = button.querySelector("[data-pick-overlay-toggle-label]");
            if (label) label.textContent = activeModeLabel;
            const adjustableSize = activeMode === "circle" || activeMode === "boxes";
            if (sizeControl) sizeControl.classList.toggle("hidden", !adjustableSize);
            if (sizeInput) sizeInput.disabled = !adjustableSize;
            if (sizeNumber) sizeNumber.disabled = !adjustableSize;
            redrawOverlays(container);
        };

        button.addEventListener("click", () => {
            const activeMode = modes.includes(button.dataset.pickOverlayMode)
                ? button.dataset.pickOverlayMode
                : defaultMode;
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
        if (colorControl) {
            colorToggle?.addEventListener("click", () => {
                setColorPanelOpen(colorToggle.getAttribute("aria-expanded") !== "true");
            });
            for (const preset of colorPresets) {
                preset.addEventListener("click", () => {
                    setColor(preset.dataset.pickOverlayColorPreset);
                });
            }
            colorDisc?.addEventListener("pointerdown", (event) => {
                event.preventDefault();
                colorDiscDragging = true;
                colorDisc.focus();
                colorDisc.setPointerCapture?.(event.pointerId);
                setColorFromDiscCoordinates(event.clientX, event.clientY);
            });
            colorDisc?.addEventListener("pointermove", (event) => {
                if (colorDiscDragging) {
                    setColorFromDiscCoordinates(event.clientX, event.clientY);
                }
            });
            colorDisc?.addEventListener("pointerup", (event) => {
                if (colorDiscDragging) {
                    setColorFromDiscCoordinates(event.clientX, event.clientY);
                }
                colorDiscDragging = false;
                if (colorDisc.hasPointerCapture?.(event.pointerId)) {
                    colorDisc.releasePointerCapture(event.pointerId);
                }
            });
            colorDisc?.addEventListener("pointercancel", () => {
                colorDiscDragging = false;
            });
            colorDisc?.addEventListener("keydown", (event) => {
                const step = event.shiftKey ? 10 : 1;
                let handled = true;
                if (event.key === "ArrowLeft") {
                    colorState.hue = (colorState.hue - step + 360) % 360;
                } else if (event.key === "ArrowRight") {
                    colorState.hue = (colorState.hue + step) % 360;
                } else if (event.key === "ArrowDown") {
                    colorState.saturation = Math.max(0, colorState.saturation - step);
                } else if (event.key === "ArrowUp") {
                    colorState.saturation = Math.min(100, colorState.saturation + step);
                } else {
                    handled = false;
                }
                if (handled) {
                    event.preventDefault();
                    commitColor(hsvToHex(
                        colorState.hue,
                        colorState.saturation,
                        colorState.brightness,
                    ));
                }
            });
            colorControl.addEventListener("keydown", (event) => {
                if (event.key === "Escape") {
                    setColorPanelOpen(false);
                    colorToggle?.focus();
                }
            });
            document.addEventListener("click", (event) => {
                if (!colorControl.contains(event.target)) setColorPanelOpen(false);
            });
            colorHex?.addEventListener("input", () => {
                colorHex.setAttribute(
                    "aria-invalid",
                    setColor(colorHex.value) ? "false" : "true",
                );
            });
            colorHex?.addEventListener("change", () => {
                if (!setColor(colorHex.value)) {
                    updateColorInterface(
                        container.dataset.pickOverlayColor || defaultPickOverlayColor,
                    );
                }
            });
            colorHex?.addEventListener("keydown", (event) => {
                if (event.key === "Enter" && setColor(colorHex.value)) {
                    event.preventDefault();
                    setColorPanelOpen(false);
                    colorToggle?.focus();
                }
            });
            for (const [channel, input] of Object.entries(colorSliders)) {
                input.addEventListener("input", () => {
                    colorState[channel] = Number(input.value);
                    commitColor(hsvToHex(
                        colorState.hue,
                        colorState.saturation,
                        colorState.brightness,
                    ));
                });
            }
            setColor(container.dataset.pickOverlayColor || defaultPickOverlayColor);
        }
        setMode(button.dataset.pickOverlayMode || defaultMode);
    };

    const initializeDisplayControls = (controls) => {
        const target = document.getElementById(controls.dataset.pickDisplayTargetId || "");
        if (!target) return;

        const displayFilter = [
            "brightness(var(--pick-level-brightness,1))",
            "contrast(var(--pick-level-contrast,1))",
        ].join(" ");
        const imageSelector = controls.dataset.pickDisplayImageSelector
            || "[data-pick-micrograph-slide] > img";
        for (const image of target.querySelectorAll(imageSelector)) {
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
        document.querySelectorAll("[data-class-average-overlay-preview]")
            .forEach(initializeClassAverageOverlay);
        document.querySelectorAll("[data-pick-overlay-toggle]").forEach(initializeOverlayToggle);
        document.querySelectorAll("[data-pick-display-controls]").forEach(initializeDisplayControls);
    };
    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", initializeAll, {once: true});
    } else {
        initializeAll();
    }
})();
