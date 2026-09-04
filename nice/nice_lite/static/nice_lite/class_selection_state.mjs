export function sortedClasses(classes, sortKey) {
    const copy = [...classes];
    const missingLast = (value) =>
        typeof value === "number" && Number.isFinite(value) ? [0, value] : [1, 0];

    copy.sort((left, right) => {
        let leftKey;
        let rightKey;
        if (sortKey === "resolution") {
            leftKey = missingLast(left.resolution);
            rightKey = missingLast(right.resolution);
        } else if (sortKey === "population") {
            const leftValue = missingLast(left.population);
            const rightValue = missingLast(right.population);
            leftKey = [leftValue[0], -leftValue[1]];
            rightKey = [rightValue[0], -rightValue[1]];
        } else {
            leftKey = [0, left.class_id];
            rightKey = [0, right.class_id];
        }

        if (leftKey[0] !== rightKey[0]) return leftKey[0] - rightKey[0];
        if (leftKey[1] !== rightKey[1]) return leftKey[1] - rightKey[1];
        return left.class_id - right.class_id;
    });
    return copy;
}

export function toggleClass(selected, classId) {
    const result = new Set(selected);
    if (result.has(classId)) result.delete(classId);
    else result.add(classId);
    return result;
}

export function applyVisualRange(selected, orderedClassIds, anchorId, targetId) {
    const anchorPosition = orderedClassIds.indexOf(anchorId);
    const targetPosition = orderedClassIds.indexOf(targetId);
    if (anchorPosition < 0 || targetPosition < 0) {
        return toggleClass(selected, targetId);
    }

    const result = new Set(selected);
    const shouldSelect = selected.has(anchorId);
    const start = Math.min(anchorPosition, targetPosition);
    const end = Math.max(anchorPosition, targetPosition);
    for (const classId of orderedClassIds.slice(start, end + 1)) {
        if (shouldSelect) result.add(classId);
        else result.delete(classId);
    }
    return result;
}

export function normalizedClassIds(values, knownIds) {
    const known = new Set(knownIds);
    if (!Array.isArray(values)) return [];
    return [...new Set(
        values.filter((value) => Number.isInteger(value) && known.has(value)),
    )].sort((left, right) => left - right);
}

export function wheelAdjustedTileSize(current, deltaY, minimum, maximum, step) {
    const value = Number(current);
    const lower = Number(minimum);
    const upper = Number(maximum);
    const increment = Math.abs(Number(step));
    if (
        ![value, deltaY, lower, upper, increment].every(Number.isFinite)
        || increment === 0
    ) {
        return value;
    }
    if (deltaY === 0) return Math.min(upper, Math.max(lower, value));
    const next = value + (deltaY < 0 ? increment : -increment);
    return Math.min(upper, Math.max(lower, next));
}
