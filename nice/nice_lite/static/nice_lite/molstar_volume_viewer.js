const DEFAULT_SURFACE_COLOR = 0x808080;
const DEFAULT_CONTOUR_FRACTION = 0.3;
const STARTING_ZOOM_FACTOR = 2;
const BACKGROUND_COLORS = {
    black: 0x000000,
    white: 0xffffff,
};

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

function initialIsovalue(option) {
    const range = volumeDensityRange(option);
    if (range) {
        return {
            type: "absolute",
            value: range.minimum
                + DEFAULT_CONTOUR_FRACTION * (range.maximum - range.minimum),
            color: DEFAULT_SURFACE_COLOR,
        };
    }
    return {
        type: "relative",
        value: 1.5,
        color: DEFAULT_SURFACE_COLOR,
    };
}

function formatDensity(value) {
    return Number(value.toPrecision(5)).toString();
}

function setIsovalueDisplay(control, value) {
    if (control instanceof HTMLInputElement) {
        control.value = value;
    } else {
        control.textContent = value;
    }
}

function setIsovalueDisplayDisabled(control, disabled) {
    if (control instanceof HTMLInputElement) {
        control.disabled = disabled;
    }
}

function configureIsovalueControl(option, input, display) {
    const isovalue = initialIsovalue(option);
    if (isovalue.type !== "absolute") {
        input.disabled = true;
        setIsovalueDisplay(display, "n/a");
        setIsovalueDisplayDisabled(display, true);
        return isovalue;
    }

    const range = volumeDensityRange(option);
    input.min = String(range.minimum);
    input.max = String(range.maximum);
    input.step = String((range.maximum - range.minimum) / 1000);
    input.value = String(isovalue.value);
    setIsovalueDisplay(display, formatDensity(isovalue.value));
    setIsovalueDisplayDisabled(display, false);
    return isovalue;
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
    const isovalueOutput = root.querySelector("[data-volume-isovalue-output]");
    const isovalueDisplay = isovalueText || isovalueOutput;
    const resetButton = root.querySelector("[data-volume-reset]");
    const viewport = root.querySelector("[data-volume-viewport]");
    const host = root.querySelector("[data-volume-molstar]");
    const metadata = root.querySelector("[data-volume-metadata]");
    const loadingIndicator = root.querySelector("[data-volume-loading]");
    if (
        !sourceSelect
        || !backgroundSelect
        || !isovalueInput
        || !isovalueDisplay
        || !resetButton
        || !viewport
        || !host
        || !metadata
    ) {
        return;
    }
    if (!window.molstar?.Viewer) {
        return;
    }

    const viewer = await window.molstar.Viewer.create(host, {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: false,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        collapseRightPanel: false,
        viewportShowExpand: false,
        viewportShowToggleFullscreen: true,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
        viewportShowTrajectoryControls: false,
        viewportBackgroundColor: "#ffffff",
        powerPreference: "high-performance",
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
        const color = BACKGROUND_COLORS[backgroundSelect.value]
            ?? BACKGROUND_COLORS.white;
        viewport.style.backgroundColor = color === BACKGROUND_COLORS.white
            ? "#ffffff"
            : "#000000";
        viewer.plugin.canvas3d?.setProps({renderer: {backgroundColor: color}});
    };

    let isovalueUpdateTimer = null;
    let isovalueUpdatePromise = Promise.resolve();
    const applyIsovalue = (value) => {
        isovalueUpdatePromise = isovalueUpdatePromise
            .catch(() => {})
            .then(() => updateMolstarIsovalue(viewer, value));
        return isovalueUpdatePromise;
    };

    const scheduleIsovalueUpdate = () => {
        const value = Number(isovalueInput.value);
        if (!Number.isFinite(value)) {
            return;
        }
        setIsovalueDisplay(isovalueDisplay, formatDensity(value));
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
        setIsovalueDisplayDisabled(isovalueDisplay, true);
        resetButton.disabled = true;
        metadata.textContent = formatVolumeMetadata(option);
        loadingIndicator?.classList.remove("hidden");
        window.clearTimeout(isovalueUpdateTimer);
        await isovalueUpdatePromise.catch(() => {});
        const isovalue = configureIsovalueControl(
            option,
            isovalueInput,
            isovalueDisplay,
        );
        try {
            await viewer.plugin.clear();
            await viewer.loadVolumeFromUrl(
                {url: option.value, format: "ccp4", isBinary: true},
                [isovalue],
                {entryId: option.dataset.volumeName || "abinitio3D"},
            );
            applyBackground();
            viewer.plugin.canvas3d?.requestCameraReset({
                durationMs: 0,
                snapshot: startingCameraSnapshot,
            });
            isovalueInput.disabled = isovalue.type !== "absolute";
            setIsovalueDisplayDisabled(isovalueDisplay, isovalue.type !== "absolute");
        } catch (error) {
            console.error("Unable to load the Batch volume in Mol*.", error);
        } finally {
            sourceSelect.disabled = false;
            resetButton.disabled = false;
            loadingIndicator?.classList.add("hidden");
        }
    };

    backgroundSelect.addEventListener("change", applyBackground);
    isovalueInput.addEventListener("input", scheduleIsovalueUpdate);
    isovalueText?.addEventListener("input", () => {
        const value = Number(isovalueText.value);
        if (!Number.isFinite(value)) {
            return;
        }
        isovalueInput.value = String(value);
        setIsovalueDisplay(isovalueDisplay, formatDensity(Number(isovalueInput.value)));
        scheduleIsovalueUpdate();
    });
    resetButton.addEventListener("click", async () => {
        const isovalue = configureIsovalueControl(
            selectedVolumeOption(sourceSelect),
            isovalueInput,
            isovalueDisplay,
        );
        if (isovalue.type === "absolute") {
            try {
                await applyIsovalue(isovalue.value);
                isovalueInput.disabled = false;
                setIsovalueDisplayDisabled(isovalueDisplay, false);
            } catch (error) {
                console.error("Unable to reset the Mol* iso value.", error);
            }
        }
        viewer.plugin.canvas3d?.requestCameraReset({
            durationMs: 250,
            snapshot: startingCameraSnapshot,
        });
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
