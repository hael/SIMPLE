const VOLUME_COLORS = {
    viridis: [0.122, 0.647, 0.529],
    plasma: [0.798, 0.280, 0.470],
    inferno: [0.735, 0.216, 0.330],
    magma: [0.716, 0.215, 0.475],
    bone: [0.690, 0.750, 0.770],
    coolwarm: [0.865, 0.394, 0.322],
};

const BACKGROUND_COLORS = {
    black: [0.067, 0.067, 0.067],
    white: [1.0, 1.0, 1.0],
};
const MIN_CAMERA_DISTANCE = 0.35;
const MAX_CAMERA_DISTANCE = 6;

const VERTEX_SHADER = `#version 300 es
in vec2 aPosition;
out vec2 vPosition;

void main() {
    vPosition = aPosition;
    gl_Position = vec4(aPosition, 0.0, 1.0);
}
`;

const FRAGMENT_SHADER = `#version 300 es
precision highp float;
precision highp sampler3D;

in vec2 vPosition;
out vec4 outputColor;

uniform sampler3D uVolume;
uniform vec3 uScale;
uniform vec3 uTextureStep;
uniform mat3 uRotation;
uniform mat3 uInverseRotation;
uniform float uThreshold;
uniform float uAspect;
uniform float uTanHalfFov;
uniform float uCameraDistance;
uniform vec3 uSurfaceColor;
uniform vec3 uBackgroundColor;
uniform float uLightBackground;

vec2 intersectBox(vec3 origin, vec3 direction, vec3 boxMinimum, vec3 boxMaximum) {
    vec3 safeDirection = vec3(
        abs(direction.x) < 0.000001 ? 0.000001 : direction.x,
        abs(direction.y) < 0.000001 ? 0.000001 : direction.y,
        abs(direction.z) < 0.000001 ? 0.000001 : direction.z
    );
    vec3 first = (boxMinimum - origin) / safeDirection;
    vec3 second = (boxMaximum - origin) / safeDirection;
    vec3 nearer = min(first, second);
    vec3 farther = max(first, second);
    return vec2(
        max(max(nearer.x, nearer.y), nearer.z),
        min(min(farther.x, farther.y), farther.z)
    );
}

float densityAt(vec3 localPosition) {
    vec3 texturePosition = clamp((localPosition / uScale) + 0.5, 0.0, 1.0);
    return texture(uVolume, texturePosition).r;
}

vec3 surfaceNormal(vec3 localPosition) {
    vec3 texturePosition = clamp((localPosition / uScale) + 0.5, 0.0, 1.0);
    float gradientX = texture(uVolume, texturePosition + vec3(uTextureStep.x, 0.0, 0.0)).r
        - texture(uVolume, texturePosition - vec3(uTextureStep.x, 0.0, 0.0)).r;
    float gradientY = texture(uVolume, texturePosition + vec3(0.0, uTextureStep.y, 0.0)).r
        - texture(uVolume, texturePosition - vec3(0.0, uTextureStep.y, 0.0)).r;
    float gradientZ = texture(uVolume, texturePosition + vec3(0.0, 0.0, uTextureStep.z)).r
        - texture(uVolume, texturePosition - vec3(0.0, 0.0, uTextureStep.z)).r;
    return normalize(vec3(gradientX, gradientY, gradientZ) / uScale);
}

void main() {
    vec3 rayOriginWorld = vec3(0.0, 0.0, uCameraDistance);
    vec3 rayDirectionWorld = normalize(vec3(
        vPosition.x * uAspect * uTanHalfFov,
        vPosition.y * uTanHalfFov,
        -1.0
    ));
    vec3 rayOrigin = uInverseRotation * rayOriginWorld;
    vec3 rayDirection = uInverseRotation * rayDirectionWorld;
    vec2 intersection = intersectBox(
        rayOrigin,
        rayDirection,
        -0.5 * uScale,
        0.5 * uScale
    );

    float rayStart = max(intersection.x, 0.0);
    float rayEnd = intersection.y;
    if (rayEnd <= rayStart) {
        outputColor = vec4(uBackgroundColor, 1.0);
        return;
    }

    float previousDistance = rayStart;
    float previousValue = densityAt(rayOrigin + rayDirection * rayStart) - uThreshold;
    bool foundSurface = false;
    float surfaceDistance = rayStart;

    for (int step = 1; step <= 256; step += 1) {
        float currentDistance = mix(rayStart, rayEnd, float(step) / 256.0);
        float currentValue = densityAt(rayOrigin + rayDirection * currentDistance) - uThreshold;
        if (
            (previousValue <= 0.0 && currentValue >= 0.0)
            || (previousValue >= 0.0 && currentValue <= 0.0)
        ) {
            float lowerDistance = previousDistance;
            float upperDistance = currentDistance;
            float lowerValue = previousValue;
            for (int refinement = 0; refinement < 5; refinement += 1) {
                float middleDistance = 0.5 * (lowerDistance + upperDistance);
                float middleValue = densityAt(
                    rayOrigin + rayDirection * middleDistance
                ) - uThreshold;
                if (
                    (lowerValue <= 0.0 && middleValue <= 0.0)
                    || (lowerValue >= 0.0 && middleValue >= 0.0)
                ) {
                    lowerDistance = middleDistance;
                    lowerValue = middleValue;
                } else {
                    upperDistance = middleDistance;
                }
            }
            surfaceDistance = 0.5 * (lowerDistance + upperDistance);
            foundSurface = true;
            break;
        }
        previousDistance = currentDistance;
        previousValue = currentValue;
    }

    if (!foundSurface) {
        outputColor = vec4(uBackgroundColor, 1.0);
        return;
    }

    vec3 localPosition = rayOrigin + rayDirection * surfaceDistance;
    vec3 normal = normalize(uRotation * surfaceNormal(localPosition));
    if (dot(normal, rayDirectionWorld) > 0.0) {
        normal = -normal;
    }
    vec3 lightDirection = normalize(vec3(0.45, 0.65, 1.0));
    float diffuse = max(dot(normal, lightDirection), 0.0);
    float facing = max(dot(normal, -rayDirectionWorld), 0.0);
    float specular = pow(max(dot(reflect(-lightDirection, normal), -rayDirectionWorld), 0.0), 24.0);
    float rim = pow(1.0 - facing, 2.0);
    vec3 shaded = uSurfaceColor * (0.24 + 0.72 * diffuse)
        + vec3(0.24 * specular)
        + uSurfaceColor * (0.12 * rim);
    shaded *= mix(1.0, 0.76, uLightBackground);
    outputColor = vec4(clamp(shaded, 0.0, 1.0), 1.0);
}
`;

function compileShader(gl, shaderType, source) {
    const shader = gl.createShader(shaderType);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        const message = gl.getShaderInfoLog(shader) || "unknown shader error";
        gl.deleteShader(shader);
        throw new Error(message);
    }
    return shader;
}

function createProgram(gl) {
    const program = gl.createProgram();
    const vertexShader = compileShader(gl, gl.VERTEX_SHADER, VERTEX_SHADER);
    const fragmentShader = compileShader(gl, gl.FRAGMENT_SHADER, FRAGMENT_SHADER);
    gl.attachShader(program, vertexShader);
    gl.attachShader(program, fragmentShader);
    gl.linkProgram(program);
    gl.deleteShader(vertexShader);
    gl.deleteShader(fragmentShader);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
        const message = gl.getProgramInfoLog(program) || "unknown WebGL link error";
        gl.deleteProgram(program);
        throw new Error(message);
    }
    return program;
}

function parseNumberList(value, expectedLength) {
    const numbers = String(value || "").split(",").map(Number);
    return numbers.length === expectedLength && numbers.every(Number.isFinite)
        ? numbers
        : null;
}

function rotationRows(yaw, pitch) {
    const cosineYaw = Math.cos(yaw);
    const sineYaw = Math.sin(yaw);
    const cosinePitch = Math.cos(pitch);
    const sinePitch = Math.sin(pitch);
    return [
        [cosineYaw, sineYaw * sinePitch, sineYaw * cosinePitch],
        [0, cosinePitch, -sinePitch],
        [-sineYaw, cosineYaw * sinePitch, cosineYaw * cosinePitch],
    ];
}

function transpose(rows) {
    return rows[0].map((_, column) => rows.map((row) => row[column]));
}

function columnMajor(rows) {
    return new Float32Array([
        rows[0][0], rows[1][0], rows[2][0],
        rows[0][1], rows[1][1], rows[2][1],
        rows[0][2], rows[1][2], rows[2][2],
    ]);
}

function transformDirection(rows, direction) {
    return rows.map((row) => (
        row[0] * direction[0]
        + row[1] * direction[1]
        + row[2] * direction[2]
    ));
}

function drawOrientationAxes(canvas, rows) {
    const context = canvas.getContext("2d");
    context.clearRect(0, 0, canvas.width, canvas.height);
    const axes = [
        {label: "x", color: "#ff5f56", vector: [1, 0, 0]},
        {label: "y", color: "#32d74b", vector: [0, 1, 0]},
        {label: "z", color: "#64a8ff", vector: [0, 0, 1]},
    ].map((axis) => ({
        ...axis,
        projected: transformDirection(rows, axis.vector),
    })).sort((left, right) => left.projected[2] - right.projected[2]);
    const origin = [42, 54];
    const length = 28;
    context.lineWidth = 2;
    context.font = "bold 11px sans-serif";
    for (const axis of axes) {
        const endpoint = [
            origin[0] + axis.projected[0] * length,
            origin[1] - axis.projected[1] * length,
        ];
        context.strokeStyle = axis.color;
        context.fillStyle = axis.color;
        context.beginPath();
        context.moveTo(...origin);
        context.lineTo(...endpoint);
        context.stroke();
        context.fillText(axis.label, endpoint[0] + 3, endpoint[1] - 3);
    }
}

function initializeVolumeViewer(viewer) {
    const canvas = viewer.querySelector("[data-volume-canvas]");
    const orientationCanvas = viewer.querySelector("[data-volume-orientation]");
    const sourceSelect = viewer.querySelector("[data-volume-source]");
    const colormapSelect = viewer.querySelector("[data-volume-colormap]");
    const backgroundSelect = viewer.querySelector("[data-volume-background]");
    const thresholdInput = viewer.querySelector("[data-volume-threshold]");
    const thresholdOutput = viewer.querySelector("[data-volume-threshold-output]");
    const metadata = viewer.querySelector("[data-volume-metadata]");
    const status = viewer.querySelector("[data-volume-status]");
    const resetButton = viewer.querySelector("[data-volume-reset]");
    if (
        !canvas
        || !orientationCanvas
        || !sourceSelect
        || !colormapSelect
        || !backgroundSelect
        || !thresholdInput
        || !thresholdOutput
        || !metadata
        || !status
        || !resetButton
    ) {
        return;
    }

    const gl = canvas.getContext("webgl2", {
        alpha: false,
        antialias: true,
        powerPreference: "high-performance",
    });
    if (!gl) {
        status.textContent = "3D volume viewing requires a browser with WebGL 2 enabled.";
        status.classList.add("text-streamerror");
        return;
    }

    let program;
    try {
        program = createProgram(gl);
    } catch (error) {
        status.textContent = `Unable to initialize the 3D renderer: ${error.message}`;
        status.classList.add("text-streamerror");
        return;
    }

    const positions = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, positions);
    gl.bufferData(
        gl.ARRAY_BUFFER,
        new Float32Array([-1, -1, 3, -1, -1, 3]),
        gl.STATIC_DRAW,
    );
    const positionLocation = gl.getAttribLocation(program, "aPosition");
    gl.enableVertexAttribArray(positionLocation);
    gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

    const uniforms = Object.fromEntries([
        "uScale",
        "uTextureStep",
        "uRotation",
        "uInverseRotation",
        "uThreshold",
        "uAspect",
        "uTanHalfFov",
        "uCameraDistance",
        "uSurfaceColor",
        "uBackgroundColor",
        "uLightBackground",
        "uVolume",
    ].map((name) => [name, gl.getUniformLocation(program, name)]));

    const state = {
        texture: null,
        dimensions: null,
        sourceDimensions: null,
        voxelSize: null,
        valueRange: null,
        yaw: -0.60,
        pitch: 0.35,
        distance: 2.2,
        pointerId: null,
        pointerPosition: null,
        requestController: null,
    };

    const updateThresholdOutput = () => {
        const percentage = Number(thresholdInput.value);
        if (state.valueRange) {
            const value = state.valueRange[0]
                + (percentage / 100) * (state.valueRange[1] - state.valueRange[0]);
            thresholdOutput.textContent = `${value.toFixed(4)} (${percentage}%)`;
            thresholdInput.setAttribute(
                "aria-valuetext",
                `${value.toFixed(4)} density units`,
            );
        } else {
            thresholdOutput.textContent = `${percentage}%`;
        }
    };

    const render = () => {
        if (!state.texture || !state.dimensions || !state.voxelSize) {
            return;
        }
        const devicePixelRatio = Math.min(window.devicePixelRatio || 1, 2);
        const displayWidth = Math.max(1, Math.round(canvas.clientWidth * devicePixelRatio));
        const displayHeight = Math.max(1, Math.round(canvas.clientHeight * devicePixelRatio));
        if (canvas.width !== displayWidth || canvas.height !== displayHeight) {
            canvas.width = displayWidth;
            canvas.height = displayHeight;
        }

        const physicalSize = state.sourceDimensions.map((dimension, index) => (
            dimension * state.voxelSize[index]
        ));
        const maximumSize = Math.max(...physicalSize);
        const scale = physicalSize.map((size) => size / maximumSize);
        const rows = rotationRows(state.yaw, state.pitch);
        const surfaceColor = VOLUME_COLORS[colormapSelect.value] || VOLUME_COLORS.viridis;
        const backgroundColor = BACKGROUND_COLORS[backgroundSelect.value]
            || BACKGROUND_COLORS.black;
        const lightBackground = backgroundSelect.value === "white";

        gl.viewport(0, 0, canvas.width, canvas.height);
        gl.useProgram(program);
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_3D, state.texture);
        gl.uniform1i(uniforms.uVolume, 0);
        gl.uniform3fv(uniforms.uScale, scale);
        gl.uniform3fv(
            uniforms.uTextureStep,
            state.dimensions.map((dimension) => 1 / dimension),
        );
        gl.uniformMatrix3fv(uniforms.uRotation, false, columnMajor(rows));
        gl.uniformMatrix3fv(
            uniforms.uInverseRotation,
            false,
            columnMajor(transpose(rows)),
        );
        gl.uniform1f(uniforms.uThreshold, Number(thresholdInput.value) / 100);
        gl.uniform1f(uniforms.uAspect, canvas.width / canvas.height);
        gl.uniform1f(uniforms.uTanHalfFov, Math.tan(Math.PI / 8));
        gl.uniform1f(uniforms.uCameraDistance, state.distance);
        gl.uniform3fv(uniforms.uSurfaceColor, surfaceColor);
        gl.uniform3fv(uniforms.uBackgroundColor, backgroundColor);
        gl.uniform1f(uniforms.uLightBackground, lightBackground ? 1 : 0);
        gl.drawArrays(gl.TRIANGLES, 0, 3);
        drawOrientationAxes(orientationCanvas, rows);
    };

    const loadSelectedVolume = async () => {
        state.requestController?.abort();
        const requestController = new AbortController();
        state.requestController = requestController;
        status.classList.remove("text-streamerror");
        status.textContent = "loading volume…";
        metadata.textContent = "";
        try {
            const response = await fetch(sourceSelect.value, {
                credentials: "same-origin",
                signal: requestController.signal,
            });
            if (!response.ok) {
                throw new Error(`server returned ${response.status}`);
            }
            const dimensions = parseNumberList(
                response.headers.get("X-Volume-Dimensions"),
                3,
            );
            const sourceDimensions = parseNumberList(
                response.headers.get("X-Volume-Source-Dimensions"),
                3,
            );
            const voxelSize = parseNumberList(
                response.headers.get("X-Volume-Voxel-Size"),
                3,
            );
            const minimum = Number(response.headers.get("X-Volume-Value-Min"));
            const maximum = Number(response.headers.get("X-Volume-Value-Max"));
            if (
                !dimensions
                || !sourceDimensions
                || !voxelSize
                || !Number.isFinite(minimum)
                || !Number.isFinite(maximum)
                || maximum < minimum
            ) {
                throw new Error("volume metadata is invalid");
            }
            const bytes = new Uint8Array(await response.arrayBuffer());
            const expectedLength = dimensions.reduce((product, dimension) => (
                product * dimension
            ), 1);
            if (bytes.length !== expectedLength) {
                throw new Error("volume texture is incomplete");
            }

            if (state.texture) {
                gl.deleteTexture(state.texture);
            }
            const texture = gl.createTexture();
            gl.bindTexture(gl.TEXTURE_3D, texture);
            gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
            gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
            gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
            gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
            gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
            gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_R, gl.CLAMP_TO_EDGE);
            gl.texImage3D(
                gl.TEXTURE_3D,
                0,
                gl.R8,
                dimensions[0],
                dimensions[1],
                dimensions[2],
                0,
                gl.RED,
                gl.UNSIGNED_BYTE,
                bytes,
            );

            state.texture = texture;
            state.dimensions = dimensions;
            state.sourceDimensions = sourceDimensions;
            state.voxelSize = voxelSize;
            state.valueRange = [minimum, maximum];
            updateThresholdOutput();
            metadata.textContent = [
                `${sourceDimensions.join(" × ")} voxels`,
                `${voxelSize.map((axis) => axis.toFixed(3)).join(" × ")} Å/voxel`,
                `texture ${dimensions.join(" × ")}`,
            ].join(" · ");
            status.textContent = `showing ${sourceSelect.selectedOptions[0]?.dataset.volumeName || "volume"}`;
            render();
        } catch (error) {
            if (error.name === "AbortError") {
                return;
            }
            status.textContent = `Unable to load the volume: ${error.message}`;
            status.classList.add("text-streamerror");
        }
    };

    thresholdInput.addEventListener("input", () => {
        updateThresholdOutput();
        render();
    });
    colormapSelect.addEventListener("change", render);
    backgroundSelect.addEventListener("change", render);
    sourceSelect.addEventListener("change", loadSelectedVolume);
    resetButton.addEventListener("click", () => {
        state.yaw = -0.60;
        state.pitch = 0.35;
        state.distance = 2.2;
        render();
        canvas.focus();
    });

    canvas.addEventListener("pointerdown", (event) => {
        state.pointerId = event.pointerId;
        state.pointerPosition = [event.clientX, event.clientY];
        canvas.setPointerCapture(event.pointerId);
    });
    canvas.addEventListener("pointermove", (event) => {
        if (state.pointerId !== event.pointerId || !state.pointerPosition) {
            return;
        }
        const deltaX = event.clientX - state.pointerPosition[0];
        const deltaY = event.clientY - state.pointerPosition[1];
        state.pointerPosition = [event.clientX, event.clientY];
        state.yaw += deltaX * 0.008;
        state.pitch = Math.max(
            -Math.PI * 0.49,
            Math.min(Math.PI * 0.49, state.pitch + deltaY * 0.008),
        );
        render();
    });
    const finishPointer = (event) => {
        if (state.pointerId !== event.pointerId) {
            return;
        }
        state.pointerId = null;
        state.pointerPosition = null;
    };
    canvas.addEventListener("pointerup", finishPointer);
    canvas.addEventListener("pointercancel", finishPointer);
    canvas.addEventListener("wheel", (event) => {
        event.preventDefault();
        state.distance = Math.max(
            MIN_CAMERA_DISTANCE,
            Math.min(MAX_CAMERA_DISTANCE, state.distance * Math.exp(event.deltaY * 0.001)),
        );
        render();
    }, {passive: false});
    canvas.addEventListener("keydown", (event) => {
        const rotationStep = 0.10;
        if (event.key === "ArrowLeft") state.yaw -= rotationStep;
        else if (event.key === "ArrowRight") state.yaw += rotationStep;
        else if (event.key === "ArrowUp") state.pitch -= rotationStep;
        else if (event.key === "ArrowDown") state.pitch += rotationStep;
        else if (event.key === "+" || event.key === "=") state.distance *= 0.9;
        else if (event.key === "-" || event.key === "_") state.distance *= 1.1;
        else return;
        event.preventDefault();
        state.pitch = Math.max(-Math.PI * 0.49, Math.min(Math.PI * 0.49, state.pitch));
        state.distance = Math.max(
            MIN_CAMERA_DISTANCE,
            Math.min(MAX_CAMERA_DISTANCE, state.distance),
        );
        render();
    });

    new ResizeObserver(render).observe(canvas);
    updateThresholdOutput();
    loadSelectedVolume();
}

document.querySelectorAll("[data-volume-viewer]").forEach(initializeVolumeViewer);
