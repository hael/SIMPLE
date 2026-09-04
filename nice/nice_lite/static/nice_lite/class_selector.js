import {
    applyVisualRange,
    normalizedClassIds,
    sortedClasses,
    toggleClass,
    wheelAdjustedTileSize,
} from "./class_selection_state.mjs";

const root = document.querySelector("[data-class-selector]");
const dataNode = document.getElementById("batch-class-selector-data");
if (root && dataNode) {
    const viewer = JSON.parse(dataNode.textContent);
    const classes = Array.isArray(viewer.classes) ? viewer.classes : [];
    const knownIds = classes.map((entry) => entry.class_id);
    const originalSelection = new Set(
        normalizedClassIds(viewer.initial_selected_class_ids, knownIds),
    );
    const tiles = new Map(
        [...root.querySelectorAll("[data-class-selector-tile]")].map((tile) => [
            Number(tile.dataset.classId),
            tile,
        ]),
    );
    const grid = root.querySelector("[data-class-selector-grid]");
    const classCount = root.querySelector("[data-selected-class-count]");
    const particleCount = root.querySelector("[data-selected-particle-count]");
    const saveStatus = root.querySelector("[data-class-selector-save-status]");
    const selectionPayload = root.querySelector("[name='selected_class_ids']");
    const sortControl = root.querySelector("[data-class-selector-sort]");
    const tileSizeControl = root.querySelector("[data-class-selector-size]");
    const brightnessControl = root.querySelector("[data-class-selector-brightness]");
    const contrastControl = root.querySelector("[data-class-selector-contrast]");

    function loadStoredSelection() {
        try {
            const rawValue = window.localStorage.getItem(viewer.selection_storage_key);
            if (rawValue === null) return new Set(originalSelection);
            return new Set(normalizedClassIds(JSON.parse(rawValue), knownIds));
        } catch {
            return new Set(originalSelection);
        }
    }

    let selected = loadStoredSelection();
    let sortKey = "class";
    let anchorId = null;
    let saveTimer = null;

    function orderedEntries() {
        return sortedClasses(classes, sortKey);
    }

    function updateSummary() {
        let selectedParticles = 0;
        let totalParticles = 0;
        for (const entry of classes) {
            const population =
                typeof entry.population === "number" && Number.isFinite(entry.population)
                    ? entry.population
                    : 0;
            totalParticles += population;
            if (selected.has(entry.class_id)) selectedParticles += population;
        }
        classCount.textContent = `${selected.size.toLocaleString()} / ${classes.length.toLocaleString()}`;
        particleCount.textContent = `${selectedParticles.toLocaleString()} / ${totalParticles.toLocaleString()}`;
    }

    function render() {
        for (const entry of orderedEntries()) {
            const tile = tiles.get(entry.class_id);
            if (!tile) continue;
            const isSelected = selected.has(entry.class_id);
            tile.setAttribute("aria-pressed", isSelected ? "true" : "false");
            tile.classList.toggle("border-streamring", isSelected);
            tile.classList.toggle("bg-streamring/10", isSelected);
            tile.classList.toggle("border-streamdivider", !isSelected);
            tile.classList.toggle("bg-streambg", !isSelected);
            tile.querySelector("[data-class-selector-check]")?.classList.toggle(
                "hidden",
                !isSelected,
            );
            grid.appendChild(tile);
        }
        updateSummary();
        selectionPayload.value = JSON.stringify(
            [...selected].sort((left, right) => left - right),
        );
    }

    function saveSelection() {
        try {
            window.localStorage.setItem(
                viewer.selection_storage_key,
                JSON.stringify([...selected].sort((left, right) => left - right)),
            );
            saveStatus.textContent = "saved in this browser";
            saveStatus.classList.remove("text-streamerror");
            saveStatus.classList.add("text-streamaccent");
        } catch {
            saveStatus.textContent = "kept for this page only";
            saveStatus.classList.remove("text-streamaccent");
            saveStatus.classList.add("text-streamerror");
        }
    }

    function scheduleSave() {
        if (saveTimer !== null) window.clearTimeout(saveTimer);
        saveStatus.textContent = "unsaved changes";
        saveStatus.classList.remove("text-streamaccent", "text-streamerror");
        saveTimer = window.setTimeout(() => {
            saveTimer = null;
            saveSelection();
        }, 350);
    }

    for (const [classId, tile] of tiles) {
        tile.addEventListener("click", (event) => {
            if (event.shiftKey && anchorId !== null) {
                selected = applyVisualRange(
                    selected,
                    orderedEntries().map((entry) => entry.class_id),
                    anchorId,
                    classId,
                );
            } else {
                selected = toggleClass(selected, classId);
                anchorId = classId;
            }
            render();
            scheduleSave();
        });
    }

    sortControl.addEventListener("change", (event) => {
        sortKey = event.target.value;
        render();
    });

    tileSizeControl.addEventListener("input", (event) => {
        root.style.setProperty("--class-tile-size", `${event.target.value}px`);
    });
    grid.addEventListener(
        "wheel",
        (event) => {
            if (!event.ctrlKey && !event.metaKey) return;
            event.preventDefault();
            const value = wheelAdjustedTileSize(
                tileSizeControl.value,
                event.deltaY,
                tileSizeControl.min,
                tileSizeControl.max,
                tileSizeControl.step,
            );
            tileSizeControl.value = String(value);
            root.style.setProperty("--class-tile-size", `${value}px`);
        },
        { passive: false },
    );
    brightnessControl.addEventListener("input", (event) => {
        root.style.setProperty("--class-tile-brightness", event.target.value / 100);
    });
    contrastControl.addEventListener("input", (event) => {
        root.style.setProperty("--class-tile-contrast", event.target.value / 100);
    });

    root.querySelector("[data-class-selector-select-all]").addEventListener("click", () => {
        selected = new Set(knownIds);
        anchorId = orderedEntries().at(-1)?.class_id ?? null;
        render();
        scheduleSave();
    });
    root.querySelector("[data-class-selector-deselect-all]").addEventListener("click", () => {
        selected = new Set();
        anchorId = null;
        render();
        scheduleSave();
    });
    root.querySelector("[data-class-selector-invert]").addEventListener("click", () => {
        selected = new Set(knownIds.filter((classId) => !selected.has(classId)));
        anchorId = null;
        render();
        scheduleSave();
    });
    root.querySelector("[data-class-selector-reset]").addEventListener("click", () => {
        selected = new Set(originalSelection);
        anchorId = null;
        render();
        scheduleSave();
    });

    render();
    saveSelection();
}
