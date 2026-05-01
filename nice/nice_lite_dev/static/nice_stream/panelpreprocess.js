let lastinteraction = Date.now();
const autoscrolls = new Map();
<<<<<<< Updated upstream
=======

const stopAutoscroll = (slider) => {
    if (autoscrolls.has(slider)) {
        clearInterval(autoscrolls.get(slider));
        autoscrolls.delete(slider);
    }
};

const startAutoscroll = (slider) => {
    const id = setInterval(() => {
        const atEnd = slider.scrollLeft + slider.clientWidth >= slider.scrollWidth - 1;
        slider.scrollLeft = atEnd ? 0 : slider.scrollLeft + slider.clientWidth;
        updateLabel(slider);
    }, 5000);
    autoscrolls.set(slider, id);
};

const updateLabel = (slider) => {
    const label = slider.closest('.flex-col')?.querySelector('[id$="_label"]');
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
>>>>>>> Stashed changes

<<<<<<< HEAD
const scrlRight = (element, event) => {
    event.preventDefault();
    const slider = element.previousElementSibling;
    const count = slider.children.length;
    if (!count) return;
    const w = slider.scrollWidth / count;
    const idx = Math.round(slider.scrollLeft / w);
    slider.scrollLeft = ((idx + 1) % count) * w;
    lastinteraction = Date.now();
<<<<<<< Updated upstream
=======
const stopAutoscroll = (slider) => {
    if (autoscrolls.has(slider)) {
        clearInterval(autoscrolls.get(slider));
        autoscrolls.delete(slider);
    }
};

const startAutoscroll = (slider) => {
    const id = setInterval(() => {
        const atEnd = slider.scrollLeft + slider.clientWidth >= slider.scrollWidth - 1;
        slider.scrollLeft = atEnd ? 0 : slider.scrollLeft + slider.clientWidth;
        updateLabel(slider);
    }, 5000);
    autoscrolls.set(slider, id);
};

const updateLabel = (slider) => {
    const label = slider.closest('.flex-col')?.querySelector('[id$="_label"]');
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

const scrlRight = (element, event) => {
    event.preventDefault();
    const slider = element.previousElementSibling;
    const count = slider.children.length;
    if (!count) return;
    const w = slider.scrollWidth / count;
    const idx = Math.round(slider.scrollLeft / w);
    slider.scrollLeft = ((idx + 1) % count) * w;
    lastinteraction = Date.now();
    stopAutoscroll(slider);
>>>>>>> a1e410fad146030f3fcbc61f288170a806ef2b04
=======
    stopAutoscroll(slider);
>>>>>>> Stashed changes
};

const scrlLeft = (element, event) => {
    event.preventDefault();
<<<<<<< Updated upstream
<<<<<<< HEAD
    document.getElementById('micrograph_slider').scrollLeft -= 200;
    lastinteraction = Date.now();
=======
    const slider = element.nextElementSibling;
    const count = slider.children.length;
    if (!count) return;
    const w = slider.scrollWidth / count;
    const idx = Math.round(slider.scrollLeft / w);
    slider.scrollLeft = ((idx - 1 + count) % count) * w;
    lastinteraction = Date.now();
    stopAutoscroll(slider);
>>>>>>> a1e410fad146030f3fcbc61f288170a806ef2b04
=======
    const slider = element.nextElementSibling;
    const count = slider.children.length;
    if (!count) return;
    const w = slider.scrollWidth / count;
    const idx = Math.round(slider.scrollLeft / w);
    slider.scrollLeft = ((idx - 1 + count) % count) * w;
    lastinteraction = Date.now();
    stopAutoscroll(slider);
>>>>>>> Stashed changes
};

const showMenu = (element, event) => {
    event.preventDefault();
    const selectmenu    = element.parentElement.parentElement.querySelector("[name='selectmenu']");
    const selectmenubox = selectmenu.querySelector("[name='selectmenubox']");
    selectmenu.style.display = 'flex';
    const selectmenurect    = selectmenu.getBoundingClientRect();
    const selectmenuboxrect = selectmenubox.getBoundingClientRect();
    selectmenubox.style.left = (event.pageX + selectmenuboxrect.width > selectmenurect.width)
        ? (selectmenurect.width  - selectmenuboxrect.width)  + 'px'
        : event.pageX + 'px';
    selectmenubox.style.top  = (event.pageY + selectmenuboxrect.height > selectmenurect.height)
        ? (selectmenurect.height - selectmenuboxrect.height) + 'px'
        : event.pageY + 'px';
    lastinteraction = Date.now();
};

const hideMenu = () => {
    for (const selectmenu of document.querySelectorAll("[name='selectmenu']")) {
        selectmenu.style.display = 'none';
    }
    lastinteraction = Date.now();
};

const updateBrightness = (element) => {
    document.querySelector(':root').style.setProperty('--preprocess-brightness', element.value / 100);
};

const updateContrast = (element) => {
    document.querySelector(':root').style.setProperty('--preprocess-contrast', element.value / 100);
};

// Build a bar histogram with a draggable cutoff line.
// labelFn converts a raw label value to the annotation display string.
const buildDraggerHistogram = (canvas, formId, fieldName, labelFn) => {
    const style  = getComputedStyle(document.body);
    const labels = JSON.parse(canvas.dataset.labels.replaceAll("'", '"'));
    const data   = JSON.parse(canvas.dataset.values.replaceAll("'", '"'));

    let savedval   = Number(canvas.dataset.savedval);
    let savedidx   = labels.indexOf(savedval);
    let dragactive = false;

    const dragger = {
        id: 'dragger',
        beforeEvent(chart, args) {
            switch (args.event.type) {
                case 'mousemove':
                    if (dragactive) {
                        const xelements = chart.getElementsAtEventForMode(
                            args.event.native, 'x', { intersect: false }, true
                        );
                        if (xelements.length > 0) {
                            const idx = xelements[0].index;
                            if (chart.config.options.plugins.annotation.annotations.cutoff.value !== idx) {
                                chart.config.options.plugins.annotation.annotations.cutoff.value         = idx;
                                chart.config.options.plugins.annotation.annotations.cutoff.label.content = labelFn(labels[idx]);
                                chart.config.options.plugins.annotation.annotations.cutoffshadow.xMin    = idx;
                                savedval = labels[idx];
                                chart.update('none');
                            }
                        }
                    }
                    break;
                case 'mouseup':
                    dragactive = false;
                    document.getElementsByName(fieldName)[0].value = savedval;
                    document.getElementById(formId).submit();
                    lastinteraction = Date.now();
                    break;
                case 'mousedown':
                    dragactive = true;
                    break;
                default:
                    break;
            }
        }
    };

    new Chart(canvas.getContext('2d'), {
        type: 'bar',
        plugins: [dragger],
        options: {
            responsive: false,
            events: ['mousedown', 'mouseup', 'mousemove', 'mouseout'],
            scales: {
                x: { display: false },
                y: {
                    display: true,
                    ticks: { font: { size: 8 }, maxTicksLimit: 3 }
                }
            },
            plugins: {
                legend: { display: false },
                annotation: {
                    annotations: {
                        cutoff: {
                            type: 'line',
                            borderWidth: 0,
                            label: {
                                display: true,
                                content: labelFn(savedval),
                                position: 'start',
                                backgroundColor: style.getPropertyValue('--color-streamline').trim(),
                                color:           style.getPropertyValue('--color-streamtext').trim(),
                                font:    { size: 8 },
                                padding: 2
                            },
                            scaleID: 'x',
                            value: savedidx,
                            z: 10,
                            enter() { canvas.style.cursor = 'grab';    return true; },
                            leave() { canvas.style.cursor = 'default'; return true; }
                        },
                        cutoffshadow: {
                            type: 'box',
                            backgroundColor: 'rgba(211, 211, 211, 0.2)',
                            borderWidth: 0,
                            xMin: savedidx,
                            xMax: labels.length - 1,
                        }
                    }
                }
            }
        },
        data: {
            labels,
            datasets: [{
                data,
<<<<<<< Updated upstream
<<<<<<< HEAD
                backgroundColor: style.getPropertyValue('--color-streamring').trim(),
=======
                backgroundColor: style.getPropertyValue('--color-streamring').trim() + '33',
>>>>>>> a1e410fad146030f3fcbc61f288170a806ef2b04
=======
                backgroundColor: style.getPropertyValue('--color-streamring').trim() + '33',
>>>>>>> Stashed changes
                borderColor:     style.getPropertyValue('--color-streambar').trim(),
                hoverOffset: 4
            }]
        }
    });
};

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

    // Initialise CSS image filters for the micrograph viewer.
    const root = document.querySelector(':root');
    root.style.setProperty('--preprocess-contrast',   1.0);
    root.style.setProperty('--preprocess-brightness', 1.0);

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

<<<<<<< Updated upstream
<<<<<<< HEAD
    // // Movies doughnut
    // const style = getComputedStyle(document.body);
    // for (const canvas of document.getElementsByClassName('movies_pie_chart')) {
    //     new Chart(canvas.getContext('2d'), {
    //         type: 'doughnut',
    //         options: {
    //             responsive: false,
    //             plugins: {
    //                 legend: {
    //                     position: 'right',
    //                     labels: { boxWidth: 10, padding: 2, font: { size: 9 } }
    //                 }
    //             }
    //         },
    //         data: {
    //             labels: ['queued', 'processed', 'rejected'],
    //             datasets: [{
    //                 data: [
    //                     Number(canvas.dataset.imported) - Number(canvas.dataset.processed) - Number(canvas.dataset.rejected),
    //                     Number(canvas.dataset.processed),
    //                     Number(canvas.dataset.rejected)
    //                 ],
    //                 backgroundColor: [
    //                     style.getPropertyValue('--color-streamring').trim(),
    //                     style.getPropertyValue('--color-streamicon').trim(),
    //                     style.getPropertyValue('--color-streamrejected').trim(),
    //                 ],
    //                 hoverOffset: 4,
    //                 borderColor: style.getPropertyValue('--color-streambar').trim()
    //             }]
    //         }
    //     });
    // }

=======
>>>>>>> a1e410fad146030f3fcbc61f288170a806ef2b04
=======
>>>>>>> Stashed changes
    // Dragger histograms
    for (const canvas of document.getElementsByClassName('ctfres_histogram')) {
        buildDraggerHistogram(canvas, 'update_preprocess_ctfres', 'ctfres', v => v + 'Å');
    }
    for (const canvas of document.getElementsByClassName('astig_histogram')) {
        buildDraggerHistogram(canvas, 'update_preprocess_astig', 'astigmatism', v => v + '%');
    }
    for (const canvas of document.getElementsByClassName('icescore_histogram')) {
        buildDraggerHistogram(canvas, 'update_preprocess_icescore', 'icescore', v => Number(v).toFixed(2));
    }

    // Fade out loading gauze
    const gauze = document.getElementById('loadinggauze');
    gauze.style.opacity = '0';
    setTimeout(() => { gauze.style.display = 'none'; }, 600);

<<<<<<< Updated upstream
<<<<<<< HEAD
=======
=======
>>>>>>> Stashed changes
    // Initialise labels and auto-advance sliders every 5 s.
    const histSlider = document.getElementById('histogram_slider');
    if (histSlider) {
        updateLabel(histSlider);
        histSlider.addEventListener('scroll', () => updateLabel(histSlider), { passive: true });
        if (histSlider.scrollWidth > histSlider.clientWidth) startAutoscroll(histSlider);
    }

    const micSlider = document.getElementById('micrograph_slider');
    if (micSlider && micSlider.scrollWidth > micSlider.clientWidth) startAutoscroll(micSlider);

<<<<<<< Updated upstream
>>>>>>> a1e410fad146030f3fcbc61f288170a806ef2b04
=======
>>>>>>> Stashed changes
}, false);

window.addEventListener('visibilitychange', () => {
    if (document.visibilityState !== 'hidden') location.reload();
});

setInterval(() => {
    if ((Date.now() - lastinteraction) > 10_000 && document.visibilityState !== 'hidden') {
        lastinteraction = Date.now();
        location.reload();
    }
}, 1000);
