let lastinteraction = Date.now();
let currentMenuContainer = null;

// Build a doughnut chart with stream* CSS colour variables.
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
                    style.getPropertyValue('--color-streamring').trim(),
                    style.getPropertyValue('--color-streamicon').trim(),
                    style.getPropertyValue('--color-streamrejected').trim()
                ],
                hoverOffset: 4
            }]
        }
    });
};

const restartProcess = (element) => {
    const confirmed = confirm("Please confirm that you wish to restart this process");
    if (confirmed) element.form.submit();
};

const stopProcess = (element) => {
    const confirmed = confirm("Please confirm that you wish to stop this process");
    if (confirmed) element.form.submit();
};

// Scroll buttons use sibling traversal: [left btn][viewport][right btn]
const scrlRight = (element, event) => {
    event.preventDefault();
    const slider = element.previousElementSibling;
    const atEnd = slider.scrollLeft + slider.clientWidth >= slider.scrollWidth - 1;
    slider.scrollLeft = atEnd ? 0 : slider.scrollLeft + 215;
    lastinteraction = Date.now();
};

const scrlLeft = (element, event) => {
    event.preventDefault();
    const slider = element.nextElementSibling;
    const atStart = slider.scrollLeft <= 0;
    slider.scrollLeft = atStart ? slider.scrollWidth - slider.clientWidth : slider.scrollLeft - 215;
    lastinteraction = Date.now();
};

const toggleTemplate = (templ) => {
    templ.classList.toggle("disabledbutton");
    const xmark = templ.querySelector('.xmark');
    xmark.classList.toggle("hidden", !templ.classList.contains("disabledbutton"));
    updateCounts();
    lastinteraction = Date.now() + 90000; // don't auto-reload for 2 minutes while user is interacting
};

const selectRefs = (element) => {
    const selected = [];
    let path = "";
    for (const pick_template of document.getElementsByClassName("picktemplate")) {
        if (!pick_template.classList.contains("disabledbutton")) {
            path = pick_template.dataset.path;
            selected.push(Number(pick_template.dataset.idx));
        }
    }
    document.getElementById("final_selection_source").value = path;
    document.getElementById("final_selection").value = selected;
    if (selected.length === 0) {
        alert("You must select at least 1 reference");
    } else {
        element.form.submit();
    }
};

const showMenu = (element, event) => {
    event.preventDefault();
    currentMenuContainer = element.closest('#picking_references_final, #picking_references_available');
    const selectmenu    = element.parentElement.parentElement.parentElement.querySelector("[name='selectmenu']");
    const selectmenubox = selectmenu.querySelector("[name='selectmenubox']");
    const sortpop     = selectmenu.querySelector("#sortpop");
    const sortres     = selectmenu.querySelector("#sortres");
    const selectabove = selectmenu.querySelector("#selectabove");
    const selectbelow = selectmenu.querySelector("#selectbelow");
    selectmenubox.style.top  = event.pageY + "px";
    selectmenubox.style.left = event.pageX + "px";
    selectmenu.style.display = "flex";
    sortpop.onclick = sortPop;
    sortres.onclick = sortRes;
    if (selectabove) selectabove.onclick = () => selectAbove(element);
    if (selectbelow) selectbelow.onclick = () => selectBelow(element);
    lastinteraction = Date.now();
};

const hideMenu = () => {
    for (const selectmenu of document.querySelectorAll("[name='selectmenu']")) {
        const sortpop     = selectmenu.querySelector("#sortpop");
        const sortres     = selectmenu.querySelector("#sortres");
        const selectabove = selectmenu.querySelector("#selectabove");
        const selectbelow = selectmenu.querySelector("#selectbelow");
        selectmenu.style.display = "none";
        sortpop.onclick = null;
        sortres.onclick = null;
        if (selectabove) selectabove.onclick = null;
        if (selectbelow) selectbelow.onclick = null;
    }
    lastinteraction = Date.now();
};

const sortPop = () => {
    if (!currentMenuContainer) return;
    const templates = Array.from(currentMenuContainer.querySelectorAll(".picktemplatecontainer"));
    templates.sort((a, b) => Number(b.dataset.pop) - Number(a.dataset.pop));
    templates.forEach(item => currentMenuContainer.appendChild(item));
    hideMenu();
};

const sortRes = () => {
    if (!currentMenuContainer) return;
    const templates = Array.from(currentMenuContainer.querySelectorAll(".picktemplatecontainer"));
    templates.sort((a, b) => Number(a.dataset.res) - Number(b.dataset.res));
    templates.forEach(item => currentMenuContainer.appendChild(item));
    hideMenu();
};

const selectAbove = (element) => {
    if (!currentMenuContainer) return;
    let threshold = true;
    for (const picktemplate of currentMenuContainer.querySelectorAll(".picktemplate")) {
        const xmark = picktemplate.querySelector('.xmark');
        picktemplate.classList.toggle("disabledbutton", !threshold);
        xmark.classList.toggle("hidden", threshold);
        if (picktemplate === element) threshold = false;
    }
    hideMenu();
    updateCounts();
};

const selectBelow = (element) => {
    if (!currentMenuContainer) return;
    let threshold = false;
    for (const picktemplate of currentMenuContainer.querySelectorAll(".picktemplate")) {
        const xmark = picktemplate.querySelector('.xmark');
        if (picktemplate === element) threshold = true;
        picktemplate.classList.toggle("disabledbutton", !threshold);
        xmark.classList.toggle("hidden", threshold);
    }
    hideMenu();
    updateCounts();
};

const updateCounts = () => {
    const clustercount  = document.querySelector("#clustercount");
    const particlecount = document.querySelector("#particlecount");
    if (clustercount == null || particlecount == null) return;
    let ncls = 0, ncls_tot = 0, nptcls = 0, nptcls_tot = 0;
    for (const pick_template of document.getElementsByClassName("picktemplate")) {
        const pop = Number(pick_template.dataset.pop);
        nptcls_tot += pop;
        ncls_tot   += 1;
        if (!pick_template.classList.contains("disabledbutton")) {
            nptcls += pop;
            ncls   += 1;
        }
    }
    clustercount.innerHTML  = ncls.toLocaleString()   + " / " + ncls_tot.toLocaleString();
    particlecount.innerHTML = nptcls.toLocaleString() + " / " + nptcls_tot.toLocaleString();
};

const box_sizes = [32, 36, 40, 48, 52, 56, 64, 66, 70, 72, 80, 84, 88, 100, 104, 108, 112, 120, 128, 130, 132,
    140, 144, 150, 160, 162, 168, 176, 180, 182, 192, 200, 208, 216, 220, 224, 240, 256, 264, 288, 300, 308, 320, 324, 336,
    338, 352, 364, 384, 400, 420, 432, 448, 450, 462, 480, 486, 500, 504, 512, 520, 528, 546, 560, 576, 588, 600, 640, 648,
    650, 660, 672, 686, 700, 702, 704, 720, 726, 728, 750, 768, 770, 784, 800, 810, 840, 882, 896, 910, 924, 936, 972, 980,
    1008, 1014, 1020, 1024, 1296, 1536, 1728, 1944, 2048, 2304, 2592, 3072, 3200, 3456, 3888, 4096, 4608, 5000, 5184, 6144,
    6250, 6400, 6912, 7776, 8192, 9216, 10240, 12288, 12500];

const updateBrightness = (element) => {
    document.querySelector(':root').style.setProperty('--genpickrefs-brightness', element.value / 100);
};

const updateContrast = (element) => {
    document.querySelector(':root').style.setProperty('--genpickrefs-contrast', element.value / 100);
};

const updateScale = (element) => {
    const scale = element.value;
    for (const el of document.querySelectorAll(".cls2D, .cls2Dimg, .picktemplate, .picktemplateimg")) {
        el.style.width  = scale + "px";
        el.style.height = scale + "px";
    }
};

const drawMask = () => {
    for (const mskcanvas of document.getElementsByClassName("mskcanvas")) {
        const mskscale = Number(mskcanvas.dataset.mskscale);
        const mskdiam  = Number(mskcanvas.dataset.mskdiam);
        const ctx = mskcanvas.getContext("2d");
        ctx.strokeStyle = getComputedStyle(document.body).getPropertyValue('--color-streamaction').trim();
        ctx.clearRect(0, 0, mskcanvas.width, mskcanvas.height);
        ctx.lineWidth = 2
        ctx.beginPath();
        ctx.arc(mskcanvas.width / 2, mskcanvas.height / 2, mskdiam * mskcanvas.width / (mskscale * 2), 0, 2 * Math.PI);
        ctx.stroke();
    }
};

window.addEventListener("load", () => {
    const root = document.querySelector(':root');
    root.style.setProperty('--genpickrefs-contrast',   1.0);
    root.style.setProperty('--genpickrefs-brightness', 1.0);

    for (const canvas of document.getElementsByClassName("particles_pie_chart")) {
        const n_imported  = Number(canvas.dataset.imported);
        const n_processed = Number(canvas.dataset.processed);
        const n_rejected  = Number(canvas.dataset.rejected);
        buildDonut(canvas, [n_imported - n_processed - n_rejected, n_processed, n_rejected], ['unclassified', 'accepted', 'rejected']);
    }

    const logtext = document.querySelector(".logtext");
    if (logtext) logtext.scrollTop = logtext.scrollHeight - logtext.clientHeight;

    drawMask();
    updateCounts();

    const gauze = document.getElementById("loadinggauze");
    gauze.style.opacity = "0";
    setTimeout(() => { gauze.style.display = "none"; }, 600);
});

window.addEventListener("visibilitychange", () => {
    if (document.visibilityState !== "hidden") location.reload();
});

// 30 s idle timeout: longer than other panels because user actively selects references.
setInterval(() => {
    if ((Date.now() - lastinteraction) > 30_000 && document.visibilityState !== "hidden") {
        lastinteraction = Date.now();
        location.reload();
    }
}, 1000);
