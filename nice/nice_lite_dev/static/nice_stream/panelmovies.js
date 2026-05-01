let lastinteraction = Date.now();
let autoscroll = null;
<<<<<<< Updated upstream
=======

const updateLabel = (slider) => {
    const label = document.getElementById('histogram_label');
    if (!label) return;
    const idx = Math.round(slider.scrollLeft / slider.clientWidth);
    const slide = slider.children[idx];
    if (!slide) return;
    const newText = slide.dataset.label ?? '';
    if (label.textContent === newText) return;
    label.style.opacity = '0';
    setTimeout(() => {
        label.textContent = newText;
        label.style.opacity = '1';
    }, 150);
};

const stopAutoscroll = () => {
    if (autoscroll !== null) {
        clearInterval(autoscroll);
        autoscroll = null;
    }
};

const scrlRight = (element, event) => {
    event.preventDefault();
    const slider = element.previousElementSibling;
    const count = slider.children.length;
    if (!count) return;
    const w = slider.scrollWidth / count;
    const idx = Math.round(slider.scrollLeft / w);
    slider.scrollLeft = ((idx + 1) % count) * w;
    lastinteraction = Date.now();
    stopAutoscroll();
};

const scrlLeft = (element, event) => {
    event.preventDefault();
    const slider = element.nextElementSibling;
    const count = slider.children.length;
    if (!count) return;
    const w = slider.scrollWidth / count;
    const idx = Math.round(slider.scrollLeft / w);
    slider.scrollLeft = ((idx - 1 + count) % count) * w;
    lastinteraction = Date.now();
    stopAutoscroll();
};

const startAutoscroll = (slider) => {
    autoscroll = setInterval(() => {
        const atEnd = slider.scrollLeft + slider.clientWidth >= slider.scrollWidth - 1;
        slider.scrollLeft = atEnd ? 0 : slider.scrollLeft + slider.clientWidth;
        updateLabel(slider);
    }, 5000);
};
>>>>>>> Stashed changes

<<<<<<< HEAD
=======
const updateLabel = (slider) => {
    const label = document.getElementById('histogram_label');
    if (!label) return;
    const idx = Math.round(slider.scrollLeft / slider.clientWidth);
    const slide = slider.children[idx];
    if (!slide) return;
    const newText = slide.dataset.label ?? '';
    if (label.textContent === newText) return;
    label.style.opacity = '0';
    setTimeout(() => {
        label.textContent = newText;
        label.style.opacity = '1';
    }, 150);
};

const stopAutoscroll = () => {
    if (autoscroll !== null) {
        clearInterval(autoscroll);
        autoscroll = null;
    }
};

const scrlRight = (element, event) => {
    event.preventDefault();
    const slider = element.previousElementSibling;
    const count = slider.children.length;
    if (!count) return;
    const w = slider.scrollWidth / count;
    const idx = Math.round(slider.scrollLeft / w);
    slider.scrollLeft = ((idx + 1) % count) * w;
    lastinteraction = Date.now();
    stopAutoscroll();
};

const scrlLeft = (element, event) => {
    event.preventDefault();
    const slider = element.nextElementSibling;
    const count = slider.children.length;
    if (!count) return;
    const w = slider.scrollWidth / count;
    const idx = Math.round(slider.scrollLeft / w);
    slider.scrollLeft = ((idx - 1 + count) % count) * w;
    lastinteraction = Date.now();
    stopAutoscroll();
};

const startAutoscroll = (slider) => {
    autoscroll = setInterval(() => {
        const atEnd = slider.scrollLeft + slider.clientWidth >= slider.scrollWidth - 1;
        slider.scrollLeft = atEnd ? 0 : slider.scrollLeft + slider.clientWidth;
        updateLabel(slider);
    }, 5000);
};

>>>>>>> a1e410fad146030f3fcbc61f288170a806ef2b04
// Build a stacked bar chart using stream* CSS colour variables.
// datasets is an array of { data, stacked } objects; the second dataset
// (if present) uses streamicon as its fill colour.
const buildBar = (canvas, labels, datasets) => {
    const style = getComputedStyle(document.body);
<<<<<<< Updated upstream
<<<<<<< HEAD
    const colorPrimary   = style.getPropertyValue('--color-streamring').trim();
=======
    const colorPrimary   = style.getPropertyValue('--color-streamring').trim() + '33';
>>>>>>> a1e410fad146030f3fcbc61f288170a806ef2b04
=======
    const colorPrimary   = style.getPropertyValue('--color-streamring').trim() + '33';
>>>>>>> Stashed changes
    const colorSecondary = style.getPropertyValue('--color-streamicon').trim();
    const colorBorder    = style.getPropertyValue('--color-streambg').trim();

    // Pad labels to 24 buckets so all charts share the same x-axis width.
    const paddedLabels = [...labels];
    while (paddedLabels.length < 24) paddedLabels.push('');

    const stacked = datasets.length > 1;

    new Chart(canvas.getContext('2d'), {
        type: 'bar',
        options: {
            responsive: false,
            maxBarThickness: 4,
            scales: {
                x: { display: false, stacked },
                y: {
                    display: true,
                    stacked,
                    ticks: { font: { size: 8 }, maxTicksLimit: 3 }
                }
            },
            plugins: { legend: { display: false } }
        },
        data: {
            labels: paddedLabels,
            datasets: datasets.map((ds, i) => ({
                data: ds.data,
                hoverOffset: 4,
                backgroundColor: i === 0 ? colorPrimary : colorSecondary,
                borderColor: colorBorder
            }))
        }
    });
};

// Parse a JSON array stored in a data-* attribute (single-quoted in template).
const parseAttr = (str) => JSON.parse(str.replaceAll("'", '"'));

window.addEventListener('load', () => {

    for (const canvas of document.getElementsByClassName('astig_histogram')) {
        buildBar(canvas, parseAttr(canvas.dataset.labels), [
            { data: parseAttr(canvas.dataset.values)  },
            { data: parseAttr(canvas.dataset.values2) }
        ]);
    }

    for (const canvas of document.getElementsByClassName('ctfres_histogram')) {
        buildBar(canvas, parseAttr(canvas.dataset.labels), [
            { data: parseAttr(canvas.dataset.values)  },
            { data: parseAttr(canvas.dataset.values2) }
        ]);
    }

    for (const canvas of document.getElementsByClassName('df_histogram')) {
        buildBar(canvas, parseAttr(canvas.dataset.labels), [
            { data: parseAttr(canvas.dataset.values)  },
            { data: parseAttr(canvas.dataset.values2) }
        ]);
    }

    for (const canvas of document.getElementsByClassName('rate_histogram')) {
        buildBar(canvas, parseAttr(canvas.dataset.labels), [
            { data: parseAttr(canvas.dataset.values) }
        ]);
    }

    // Fade out the loading gauze once all charts are built.
    const gauze = document.getElementById('loadinggauze');
    gauze.style.opacity = '0';
    setTimeout(() => { gauze.style.display = 'none'; }, 600);

<<<<<<< Updated upstream
<<<<<<< HEAD
=======
=======
>>>>>>> Stashed changes
    // Initialise label and auto-advance the timeplots slider every 5 s.
    const timeplotSlider = document.getElementById('timeplots_slider');
    if (timeplotSlider) {
        updateLabel(timeplotSlider);
        timeplotSlider.addEventListener('scroll', () => updateLabel(timeplotSlider), { passive: true });
        if (timeplotSlider.scrollWidth > timeplotSlider.clientWidth) startAutoscroll(timeplotSlider);
    }

<<<<<<< Updated upstream
>>>>>>> a1e410fad146030f3fcbc61f288170a806ef2b04
=======
>>>>>>> Stashed changes
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
