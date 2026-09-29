(() => {
  "use strict";

  const container = document.getElementById("workspace_flow_graph");
  const cardsLayer = document.getElementById("workspace_flow_cards");
  const dataElement = document.getElementById("workspace_flow_data");
  const message = document.getElementById("workspace_flow_message");
  if (!container || !dataElement) return;

  const showMessage = (text) => {
    if (!message) return;
    message.textContent = text;
    message.classList.remove("hidden");
  };

  let graph;
  try {
    graph = JSON.parse(dataElement.textContent);
  } catch (_error) {
    showMessage("flow data could not be loaded");
    return;
  }

  if (!window.cytoscape) {
    showMessage("flow view could not be loaded");
    return;
  }
  if (!Array.isArray(graph.nodes) || graph.nodes.length === 0) return;

  const theme = getComputedStyle(document.documentElement);
  const themeColor = (name, fallback) => (
    theme.getPropertyValue(`--color-${name}`).trim() || fallback
  );
  const colors = {
    line: themeColor("streamlabel", "#a3a3a3"),
  };

  const cy = window.cytoscape({
    container,
    elements: [...graph.nodes, ...(graph.edges || [])],
    layout: { name: "preset" },
    minZoom: 0.2,
    maxZoom: 2.5,
    wheelSensitivity: 1,
    boxSelectionEnabled: false,
    autoungrabify: true,
    style: [
      {
        selector: "node",
        style: {
          "width": 260,
          "height": 68,
          "shape": "round-rectangle",
          "background-opacity": 0,
          "border-width": 0,
          "label": "",
          "overlay-opacity": 0,
        },
      },
      {
        selector: "edge",
        style: {
          "width": 1.5,
          "curve-style": "taxi",
          "taxi-direction": "downward",
          "taxi-turn": "50%",
          "taxi-turn-min-distance": 12,
          "source-endpoint": "outside-to-node",
          "target-endpoint": "outside-to-node",
          "line-color": colors.line,
          "target-arrow-color": colors.line,
          "target-arrow-shape": "triangle",
          "arrow-scale": 0.8,
          "opacity": 0.75,
        },
      },
    ],
  });

  // Unhide synchronously to measure the Tailwind-rendered cards, then hide the
  // layer again before the browser paints their initial top-left positions.
  cardsLayer?.classList.remove("hidden");
  const cards = new Map();
  cardsLayer?.querySelectorAll("[data-flow-node-id]").forEach((card) => {
    cards.set(card.dataset.flowNodeId, card);
    const node = cy.$id(card.dataset.flowNodeId);
    if (node.nonempty()) {
      node.style({
        width: card.offsetWidth,
        height: card.offsetHeight,
      });
    }
  });
  cardsLayer?.classList.add("hidden");

  const workspaceId = graph.summary && graph.summary.workspaceId;
  const viewportKey = `nice_lite.workspace_flow.viewport.${workspaceId || 0}`;
  let viewportSavePending = false;
  let lastInteraction = Date.now();
  window.noteJobCardInteraction = () => {
    lastInteraction = Date.now();
  };

  const saveViewport = () => {
    try {
      sessionStorage.setItem(viewportKey, JSON.stringify({
        zoom: cy.zoom(),
        pan: cy.pan(),
      }));
    } catch (_error) {
      // Storage is optional; the graph remains usable when it is unavailable.
    }
  };

  const scheduleViewportSave = () => {
    lastInteraction = Date.now();
    if (viewportSavePending) return;
    viewportSavePending = true;
    requestAnimationFrame(() => {
      viewportSavePending = false;
      saveViewport();
    });
  };

  const restoreViewport = () => {
    try {
      const stored = JSON.parse(sessionStorage.getItem(viewportKey));
      if (
        stored
        && Number.isFinite(stored.zoom)
        && stored.pan
        && Number.isFinite(stored.pan.x)
        && Number.isFinite(stored.pan.y)
      ) {
        cy.zoom(stored.zoom);
        cy.pan(stored.pan);
        return true;
      }
    } catch (_error) {
      // Ignore malformed or unavailable session storage.
    }
    return false;
  };

  const layout = cy.layout({
    name: "dagre",
    rankDir: "TB",
    nodeSep: 36,
    edgeSep: 18,
    rankSep: 40,
    padding: 64,
    fit: false,
    animate: false,
    acyclicer: "greedy",
    sort: (first, second) => (
      Number(first.data("order") || 0) - Number(second.data("order") || 0)
    ),
  });
  let cardPositionPending = false;
  const positionCards = () => {
    cardPositionPending = false;
    const zoom = cy.zoom();
    cards.forEach((card, nodeId) => {
      const node = cy.$id(nodeId);
      if (node.empty()) return;
      const position = node.renderedPosition();
      card.style.left = `${position.x}px`;
      card.style.top = `${position.y}px`;
      card.style.transform = `translate(-50%, -50%) scale(${zoom})`;
    });
    cardsLayer?.classList.remove("hidden");
  };
  const scheduleCardPosition = () => {
    if (cardPositionPending) return;
    cardPositionPending = true;
    requestAnimationFrame(positionCards);
  };
  cy.on("render", scheduleCardPosition);
  layout.on("layoutstop", () => {
    if (!restoreViewport()) cy.fit(cy.elements(), 56);
    positionCards();
    cy.on("pan zoom", scheduleViewportSave);
  });
  layout.run();

  const zoomAroundCenter = (factor) => {
    const nextZoom = Math.max(cy.minZoom(), Math.min(cy.maxZoom(), cy.zoom() * factor));
    cy.zoom({
      level: nextZoom,
      renderedPosition: {
        x: container.clientWidth / 2,
        y: container.clientHeight / 2,
      },
    });
    scheduleViewportSave();
  };

  document.getElementById("workspace_flow_zoom_out")?.addEventListener("click", () => {
    zoomAroundCenter(0.8);
  });
  document.getElementById("workspace_flow_zoom_in")?.addEventListener("click", () => {
    zoomAroundCenter(1.25);
  });
  document.getElementById("workspace_flow_fit")?.addEventListener("click", () => {
    cy.fit(cy.elements(), 56);
    scheduleViewportSave();
  });

  if (window.ResizeObserver) {
    const resizeObserver = new ResizeObserver(() => cy.resize());
    resizeObserver.observe(container);
  } else {
    window.addEventListener("resize", () => cy.resize());
  }

  window.addEventListener("pagehide", saveViewport);
  document.addEventListener("visibilitychange", () => {
    if (document.visibilityState !== "hidden") {
      saveViewport();
      location.reload();
    }
  });

  setInterval(() => {
    if ((Date.now() - lastInteraction) > 10000 && document.visibilityState !== "hidden") {
      lastInteraction = Date.now();
      saveViewport();
      location.reload();
    }
  }, 1000);
})();
