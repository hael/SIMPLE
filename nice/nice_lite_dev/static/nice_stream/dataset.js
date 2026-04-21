let lastinteraction = Date.now();

const submitUpdate = (element, event) => {
    if (event.key === 'Enter') element.form.submit();
};

const enableDatasetRename = () => {
    const el = document.querySelector('#datasetrename');
    el.disabled = false;
    el.focus();
    lastinteraction = Date.now();
};

const enableDatasetDescription = () => {
    const el = document.querySelector('#datasetdescription');
    el.disabled = false;
    el.focus();
    lastinteraction = Date.now();
};

const enableStreamDescription = (button) => {
    const input = button.parentElement.querySelector('input[name="new_stream_description"]');
    input.disabled = false;
    input.focus();
    lastinteraction = Date.now();
};

// Read job id from data-job-id attribute; submit the enclosing form on confirm.
const deleteStream = (button) => {
    const jobid = button.dataset.jobId;
    if (confirm('Please confirm that you wish to delete stream ' + jobid)) {
        button.closest('form').submit();
    }
};

const terminateStream = (button) => {
    const jobid = button.dataset.jobId;
    if (confirm('Please confirm that you wish to terminate all processes in stream ' + jobid)) {
        button.closest('form').submit();
    }
};

// Build a doughnut chart.  Chart.js legend is suppressed — the HTML card
// legend (imported/processed/rejected swatches) serves that purpose instead.
const buildDonut = (canvas, dataValues, labels) => {
    const style = getComputedStyle(document.body);
    new Chart(canvas.getContext('2d'), {
        type: 'doughnut',
        options: {
            responsive: false,
            plugins: { legend: { display: false } }
        },
        data: {
            labels,
            datasets: [{
                data: dataValues,
                backgroundColor: [
                    style.getPropertyValue('--color-streamring').trim(),     // imported / accepted
                    style.getPropertyValue('--color-streamicon').trim(),     // processed
                    style.getPropertyValue('--color-streamrejected').trim()  // rejected
                ],
                hoverOffset: 4
            }]
        }
    });
};

window.addEventListener('load', () => {

    for (const canvas of document.getElementsByClassName('movies_pie_chart')) {
        // Skip when preprocessing stats are not yet populated (empty string = no data).
        if (!canvas.dataset.imported) continue;
        buildDonut(
            canvas,
            [
                Number(canvas.dataset.imported),
                Number(canvas.dataset.processed || 0),
                Number(canvas.dataset.rejected  || 0)
            ],
            ['imported', 'processed', 'rejected']
        );
    }

    for (const canvas of document.getElementsByClassName('particles_pie_chart')) {
        if (!canvas.dataset.imported) continue;
        const imported = Number(canvas.dataset.imported);
        const accepted = Number(canvas.dataset.accepted || 0);
        const rejected = Number(canvas.dataset.rejected || 0);
        buildDonut(
            canvas,
            [imported - accepted - rejected, accepted, rejected],
            ['unclassified', 'accepted', 'rejected']
        );
    }

}, false);

// Reload when the tab becomes visible again (data may have changed while hidden).
window.addEventListener('visibilitychange', () => {
    if (document.visibilityState !== 'hidden') location.reload();
});

// Also poll every 10 s while the user is idle.
setInterval(() => {
    if ((Date.now() - lastinteraction) > 10_000 && document.visibilityState !== 'hidden') {
        lastinteraction = Date.now();
        location.reload();
    }
}, 1000);
