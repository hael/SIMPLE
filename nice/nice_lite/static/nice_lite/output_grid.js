import { wheelAdjustedTileSize } from "./class_selection_state.mjs";

function syncTileSizeControl(control, value) {
    if (!Number.isFinite(value)) return;
    control.value = String(value);
    control.setAttribute("aria-valuetext", `${control.value} pixels`);
}

document.querySelectorAll("[data-output-grid-size]").forEach((control) => {
    const target = document.getElementById(control.dataset.outputGridTargetId || "");
    if (!target) return;

    const grids = [
        ...(target.matches("[data-output-grid]") ? [target] : []),
        ...target.querySelectorAll("[data-output-grid]"),
    ].filter((grid) => grid.dataset.outputGridSizeProperty);
    if (!grids.length) return;

    const initialGrid = grids[0];
    const tileDimension = initialGrid.dataset.outputGridTileDimension === "width"
        ? "width"
        : "height";
    const initialTile = initialGrid.querySelector("[data-output-grid-tile]:not(.hidden)")
        || initialGrid.querySelector("[data-output-grid-tile]");
    const initialSize = initialTile?.getBoundingClientRect()[tileDimension];
    if (initialSize > 0) syncTileSizeControl(control, Math.round(initialSize));

    const setTileSize = (value) => {
        const boundedValue = wheelAdjustedTileSize(
            value,
            0,
            control.min,
            control.max,
            control.step,
        );
        if (!Number.isFinite(boundedValue)) return;
        syncTileSizeControl(control, boundedValue);
        grids.forEach((grid) => {
            grid.style.setProperty(
                grid.dataset.outputGridSizeProperty,
                `${control.value}px`,
            );
        });
    };

    control.addEventListener("input", (event) => {
        setTileSize(event.target.value);
    });
    grids.forEach((grid) => {
        grid.addEventListener(
            "wheel",
            (event) => {
                if (!event.ctrlKey && !event.metaKey) return;
                event.preventDefault();
                setTileSize(wheelAdjustedTileSize(
                    control.value,
                    event.deltaY,
                    control.min,
                    control.max,
                    control.step,
                ));
            },
            { passive: false },
        );
    });
});
