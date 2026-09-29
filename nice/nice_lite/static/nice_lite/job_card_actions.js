(() => {
  "use strict";

  const jobId = (button) => button?.dataset?.jobId || "";
  const submitAfterConfirmation = (button, message) => {
    if (window.confirm(message)) button.closest("form")?.submit();
    window.noteJobCardInteraction?.();
  };

  window.deleteStream = (button) => {
    const id = jobId(button);
    submitAfterConfirmation(button, `Please confirm that you wish to delete stream ${id}`);
  };

  window.terminateStream = (button) => {
    const id = jobId(button);
    submitAfterConfirmation(
      button,
      `Please confirm that you wish to terminate all processes in stream ${id}`,
    );
  };

  window.restartStream = (button) => {
    const id = jobId(button);
    submitAfterConfirmation(button, `Please confirm that you wish to restart stream ${id}`);
  };

  window.killDeleteStream = (button) => {
    const id = jobId(button);
    submitAfterConfirmation(
      button,
      `Please confirm that you wish to kill and delete stream ${id}`,
    );
  };

  window.stopBatchJob = (button) => {
    const id = jobId(button);
    submitAfterConfirmation(button, `Please confirm that you wish to stop batch job ${id}`);
  };

  window.markBatchJobFinished = (button) => {
    const id = jobId(button);
    submitAfterConfirmation(
      button,
      `Mark batch job ${id} as finished? This changes only its recorded status.`,
    );
  };

  window.deleteBatchJob = (button) => {
    const id = jobId(button);
    submitAfterConfirmation(
      button,
      `Please confirm that you wish to permanently delete batch job ${id}. This cannot be undone.`,
    );
  };

  // Drag source for the per-job data-type icon row. Output paths are attached
  // only to icons backed by an existing finished-job artifact.
  window.dragArtifact = (event, jobIdValue, dataType) => {
    const path = event.currentTarget.dataset.artifactPath || "";
    const payload = JSON.stringify({ jobId: jobIdValue, dataType, path });
    event.dataTransfer.setData("application/json", payload);
    event.dataTransfer.setData("text/plain", payload);
    event.dataTransfer.effectAllowed = "copy";

    const svgSource = event.currentTarget.querySelector("svg");
    if (svgSource) {
      const ghost = svgSource.cloneNode(true);
      ghost.removeAttribute("class");
      ghost.style.width = "28px";
      ghost.style.height = "28px";
      ghost.style.color = getComputedStyle(event.currentTarget).color;
      ghost.style.position = "fixed";
      ghost.style.top = "-200px";
      ghost.style.left = "-200px";
      ghost.style.pointerEvents = "none";
      document.body.appendChild(ghost);
      event.dataTransfer.setDragImage(ghost, 14, 14);
      setTimeout(() => ghost.remove(), 0);
    }

    window.noteJobCardInteraction?.();
  };

  const isJobBuilderActive = () => {
    try {
      const iframe = window.parent?.document?.getElementById("job_builder_iframe");
      if (!iframe || iframe.classList.contains("hidden")) return false;
      return window.parent.getComputedStyle(iframe).display !== "none";
    } catch (_error) {
      return false;
    }
  };

  const artifactRows = document.querySelectorAll(".artifact-drag-row");
  const batchProjectCards = document.querySelectorAll("[data-batch-project-path]");
  if (artifactRows.length === 0 && batchProjectCards.length === 0) return;

  const syncArtifactRows = () => {
    const active = isJobBuilderActive();
    artifactRows.forEach((row) => {
      row.classList.toggle("hidden", !active);
    });

    batchProjectCards.forEach((card) => {
      card.draggable = active;
      card.classList.toggle("cursor-pointer", !active);
      card.classList.toggle("cursor-grab", active);
      card.classList.toggle("active:cursor-grabbing", active);
      if (active) card.setAttribute("title", "drag to input project");
      else card.removeAttribute("title");
    });
  };

  syncArtifactRows();
  setInterval(syncArtifactRows, 500);
})();
