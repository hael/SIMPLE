const DEFAULT_SURFACE_COLOR = 0x808080;
const DEFAULT_VISIBLE_VOXEL_FRACTION = 0.01;
const STARTING_ZOOM_FACTOR = 2;

function selectedVolumeOption(sourceSelect) {
    return sourceSelect.options[sourceSelect.selectedIndex] || null;
}

function datasetNumber(option, name) {
    const value = Number(option?.dataset[name]);
    return Number.isFinite(value) ? value : null;
}

function volumeDensityRange(option) {
    const minimum = datasetNumber(option, "volumeMinimum");
    const maximum = datasetNumber(option, "volumeMaximum");
    return minimum !== null && maximum !== null && maximum > minimum
        ? {minimum, maximum}
        : null;
}

function densityToInt8(option, value) {
    const range = volumeDensityRange(option);
    // Match the server's linear conversion while keeping controls in source units.
    return range
        ? Math.max(-128, Math.min(127,
            -128 + 255 * (value - range.minimum) / (range.maximum - range.minimum)))
        : value;
}

function int8ToDensity(option, value) {
    const range = volumeDensityRange(option);
    return range
        ? range.minimum + (value + 128) * (range.maximum - range.minimum) / 255
        : value;
}

function topOnePercentThreshold(values) {
    // The volume endpoint sends signed 8-bit voxels, so a histogram avoids
    // sorting or copying the full grid in the browser.
    if (!(values instanceof Int8Array) || !values.length) {
        throw new Error("The loaded volume has no 8-bit voxel data.");
    }
    const counts = new Uint32Array(256);
    for (const value of values) {
        counts[value + 128] += 1;
    }

    let remaining = Math.max(1, Math.ceil(values.length * DEFAULT_VISIBLE_VOXEL_FRACTION));
    for (let value = 127; value >= -128; value--) {
        remaining -= counts[value + 128];
        if (remaining <= 0) {
            return value;
        }
    }
}

function formatDensity(value) {
    return Number(value.toPrecision(5)).toString();
}

function configureIsovalueControl(option, input, text, value) {
    const range = volumeDensityRange(option);
    for (const control of [input, text]) {
        control.min = String(range?.minimum ?? -128);
        control.max = String(range?.maximum ?? 127);
        control.step = "any";
    }
    input.value = String(value);
    text.value = formatDensity(value);
}

async function updateMolstarIsovalue(viewer, absoluteValue) {
    const representation = viewer.plugin.managers.volume.hierarchy.current
        .volumes[0]?.representations[0]?.cell;
    const transform = representation?.transform;
    if (!transform?.params || transform.params.type?.name !== "isosurface") {
        return false;
    }

    await viewer.plugin.build().to(transform.ref).update({
        ...transform.params,
        type: {
            ...transform.params.type,
            params: {
                ...transform.params.type.params,
                isoValue: {kind: "absolute", absoluteValue},
            },
        },
    }).commit();
    return true;
}

async function addMolstarSurface(viewer, volume, absoluteValue) {
    const {registry, themes} = viewer.plugin.representation.volume;
    const isosurface = registry.get("isosurface");
    const color = themes.colorThemeRegistry.get("uniform");
    const size = themes.sizeThemeRegistry.get("uniform");
    await viewer.plugin.build().to(volume.cell).apply(
        window.molstar.lib.plugin.StateTransforms.Representation.VolumeRepresentation3D,
        {
            type: {name: "isosurface", params: {
                ...isosurface.defaultValues,
                isoValue: {kind: "absolute", absoluteValue},
            }},
            colorTheme: {name: "uniform", params: {
                ...color.defaultValues,
                value: DEFAULT_SURFACE_COLOR,
            }},
            sizeTheme: {name: "uniform", params: size.defaultValues},
        },
    ).commit();
}

function formatVolumeMetadata(option) {
    const dimensions = ["volumeWidth", "volumeHeight", "volumeDepth"]
        .map((name) => datasetNumber(option, name));
    const voxelSize = ["volumeVoxelX", "volumeVoxelY", "volumeVoxelZ"]
        .map((name) => datasetNumber(option, name));
    const details = [];
    if (dimensions.every((value) => value !== null)) {
        details.push(`${dimensions.join(" × ")} voxels`);
    }
    if (voxelSize.every((value) => value !== null)) {
        details.push(`${voxelSize.map((value) => value.toPrecision(4)).join(" × ")} Å/voxel`);
    }
    return details.join(" · ");
}

function startingCameraSnapshot(scene, camera) {
    const {center, radius} = scene.boundingSphereVisible;
    const snapshot = camera.getFocus(center, radius);
    const {position, target} = snapshot;
    if (!position || !target) {
        return snapshot;
    }
    return {
        ...snapshot,
        position: position.map((coordinate, index) => (
            target[index]
            + (coordinate - target[index]) / STARTING_ZOOM_FACTOR
        )),
    };
}

async function initializeMolstarVolumeViewer(root) {
    const sourceSelect = root.querySelector("[data-volume-source]");
    const backgroundSelect = root.querySelector("[data-volume-background]");
    const isovalueInput = root.querySelector("[data-volume-isovalue]");
    const isovalueText = root.querySelector("[data-volume-isovalue-text]");
    const resetButton = root.querySelector("[data-volume-reset]");
    const viewport = root.querySelector("[data-volume-viewport]");
    const host = root.querySelector("[data-volume-molstar]");
    const metadata = root.querySelector("[data-volume-metadata]");
    const loadingIndicator = root.querySelector("[data-volume-loading]");
    if (!window.molstar?.Viewer || [
        sourceSelect, backgroundSelect, isovalueInput, isovalueText,
        resetButton, viewport, host, metadata,
    ].some((element) => !element)) {
        return;
    }

    const viewer = await window.molstar.Viewer.create(host, {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: false,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: false,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
        viewportShowTrajectoryControls: false,
        viewportBackgroundColor: "#ffffff",
    });

    viewer.plugin.canvas3d?.setProps({
        trackball: {
            minDistance: 0.01,
            autoAdjustMinMaxDistance: {
                name: "on",
                params: {
                    minDistanceFactor: 0,
                    minDistancePadding: 0.1,
                    maxDistanceFactor: 10,
                    maxDistanceMin: 20,
                },
            },
        },
    });

    const applyBackground = () => {
        const black = backgroundSelect.value === "black";
        const color = black ? 0x000000 : 0xffffff;
        viewport.style.backgroundColor = black ? "#000000" : "#ffffff";
        viewer.plugin.canvas3d?.setProps({renderer: {backgroundColor: color}});
    };
    const resetCamera = (durationMs) => viewer.plugin.canvas3d?.requestCameraReset({
        durationMs, snapshot: startingCameraSnapshot,
    });

    let isovalueUpdateTimer = null;
    let isovalueUpdatePromise = Promise.resolve();
    let defaultDensityValue = null;
    const applyIsovalue = (value) => {
        const int8Value = densityToInt8(selectedVolumeOption(sourceSelect), value);
        isovalueUpdatePromise = isovalueUpdatePromise
            .catch(() => {})
            .then(() => updateMolstarIsovalue(viewer, int8Value));
        return isovalueUpdatePromise;
    };

    const scheduleIsovalueUpdate = (value) => {
        if (!Number.isFinite(value)) {
            return;
        }
        window.clearTimeout(isovalueUpdateTimer);
        isovalueUpdateTimer = window.setTimeout(() => {
            applyIsovalue(value).catch((error) => {
                console.error("Unable to update the Mol* iso value.", error);
            });
        }, 100);
    };

    const loadSelectedVolume = async () => {
        const option = selectedVolumeOption(sourceSelect);
        if (!option?.value) {
            return;
        }

        sourceSelect.disabled = true;
        isovalueInput.disabled = true;
        isovalueText.disabled = true;
        resetButton.disabled = true;
        metadata.textContent = formatVolumeMetadata(option);
        loadingIndicator?.classList.remove("hidden");
        window.clearTimeout(isovalueUpdateTimer);
        await isovalueUpdatePromise.catch(() => {});
        defaultDensityValue = null;
        try {
            await viewer.plugin.clear();
            await viewer.loadVolumeFromUrl(
                {url: option.value, format: "ccp4", isBinary: true},
                [],
                {entryId: option.dataset.volumeName || "solve3D"},
            );
            const volume = viewer.plugin.managers.volume.hierarchy.current.volumes[0];
            const voxels = volume?.cell?.obj?.data?.grid?.cells?.data;
            const value = int8ToDensity(
                option, topOnePercentThreshold(voxels),
            );
            configureIsovalueControl(option, isovalueInput, isovalueText, value);
            await addMolstarSurface(viewer, volume, densityToInt8(option, value));
            applyBackground();
            resetCamera(0);
            isovalueInput.disabled = false;
            isovalueText.disabled = false;
            defaultDensityValue = value;
        } catch (error) {
            console.error("Unable to load the Batch volume in Mol*.", error);
        } finally {
            sourceSelect.disabled = false;
            resetButton.disabled = defaultDensityValue === null;
            loadingIndicator?.classList.add("hidden");
        }
    };

    backgroundSelect.addEventListener("change", applyBackground);
    isovalueInput.addEventListener("input", () => {
        const value = Number(isovalueInput.value);
        isovalueText.value = formatDensity(value);
        scheduleIsovalueUpdate(value);
    });
    isovalueText.addEventListener("input", () => {
        const value = isovalueText.valueAsNumber;
        if (!Number.isFinite(value)) {
            return;
        }
        isovalueInput.value = String(value);
        scheduleIsovalueUpdate(Number(isovalueInput.value));
    });
    resetButton.addEventListener("click", async () => {
        if (defaultDensityValue === null) {
            return;
        }
        window.clearTimeout(isovalueUpdateTimer);
        configureIsovalueControl(
            selectedVolumeOption(sourceSelect), isovalueInput, isovalueText,
            defaultDensityValue,
        );
        try {
            await applyIsovalue(defaultDensityValue);
            isovalueInput.disabled = false;
            isovalueText.disabled = false;
        } catch (error) {
            console.error("Unable to reset the Mol* iso value.", error);
        }
        resetCamera(250);
    });
    sourceSelect.addEventListener("change", loadSelectedVolume);

    const resizeObserver = typeof ResizeObserver === "undefined"
        ? null
        : new ResizeObserver(() => viewer.handleResize());
    resizeObserver?.observe(host);
    window.addEventListener("pagehide", () => {
        window.clearTimeout(isovalueUpdateTimer);
        resizeObserver?.disconnect();
        viewer.dispose();
    }, {once: true});

    applyBackground();
    await loadSelectedVolume();
}

for (const root of document.querySelectorAll("[data-volume-viewer]")) {
    initializeMolstarVolumeViewer(root).catch((error) => {
        console.error("Unable to initialize the Batch Mol* viewer.", error);
    });
}
